#!/usr/bin/python

# -*- coding: utf-8 -*-
#
# This package is part of the SemNExT project at Rensselaer Polytechnic Institute.
# Copyright (c) 2015 Evan W. Patton
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Pulls in disease_filtered.csv; performs an initial ETL using Setlr and passes 
results to the UMLS endpoint to disambiguate Entrez IDs for genes.
Outputs a final disease.ttl file which contains diseases, their associated
genes, and the z-score and confidence score of the gene-disease interaction.
Run this while sourcing from the Setlr environment.
'''

author = 'Spencer Norris'

#Manage path, imports
import sys
import os
import csv
sys.path.append('./semnext')

from rdflib import Graph, ConjunctiveGraph, URIRef, RDF, Namespace
from rdflib.plugins import sparql 

from semnext.datasource.umls import UMLSDatabase
from semnext.datasource.stringdb import StringDatabase
from semnext.datasource.redrugs import RedrugsEndpoint
from semnext.datasource.ensembl import EnsemblDatabase
from semnext.datasource.bio2rdf import Bio2RDFEndpoint
from semnext.datasource.cortecon import CorteconSource
from semnext.datamodel import Gene

from SPARQLWrapper.SPARQLExceptions import EndPointInternalError
from biomart import BiomartServer
from rpy2.robjects import r as R
from rpy2.robjects.vectors import StrVector
import re

#Data source instantiations
umls_opts = {'host': '128.113.106.11', 'port': 3306, 'user': 'spencer_dev', 'password': 'nXf6-J4ma.xE3pI', 'db': 'umls2014ab'}
redrugs_opts = {'endpoint': 'http://redrugs.tw.rpi.edu/bigdata/sparql'}
string_db_opts = {'user': 'spencer_dev', 'password': 'nXf6-J4ma.xE3pI', 'db': 'string'}
#ensembl_opts = {'user': 'root', 'password': 's3cret', 'db': 'ensembl', 'host': '0.0.0.0', 'port': 32769}
bio2rdf_opts = {'endpoint': 'http://pubmed.bio2rdf.org/sparql'}
cortecon_opts = {'endpoint': 'http://localhost:8890/sparql'}
umls_db = UMLSDatabase(**umls_opts)
redrugs = RedrugsEndpoint(**redrugs_opts)
string_db = StringDatabase(**string_db_opts)
#ensembl = EnsemblDatabase(**ensembl_opts)
bio2rdf = Bio2RDFEndpoint(**bio2rdf_opts)
cortecon = CorteconSource(**cortecon_opts)

#biomart = BiomartServer( "http://www.ensembl.org/biomart" )
#humans = biomart.datasets['hsapiens_gene_ensembl']
#humans.show_filters()
#humans.show_attributes()


#Local resources
gene2ensembl_tsv = csv.reader(open('./gene2ensembl', 'r').read(), delimiter='\t')
gene2ensembl = {}
gene2ensembl_headers = ['#tax_id', 'GeneID',  'Symbol', 'LocusTag', 'Synonyms', 'dbXrefs', 'chromosome',
						'map_location', 'description', 'type_of_gene', 'Symbol_from_nomenclature_authority',
						'Full_name_from_nomenclature_authority', 'Nomenclature_status', 'Other_designations',
						'Modification_date' ]



#Misc.
TMP_PREFIX = Namespace('https://semnext.tw.rpi.edu/temporary/prefix/')
os.environ['SETLR_HOME']='/home/spencer/SemNExT/src/main/python/setlr'


'''
The different disease URIs use TMP_PREFIX:AssocatedGene
in order to tie a gene to a disease expression.
Must replace each of these.
'''
def replace_tmp_associated_genes(graph):
	#Extract all gene URIs from the graph, retrieve real URI for each
	tmp_gene_uris = set(str(gene) for gene, p, o in graph.triples( (None, RDF.type, URIRef("https://semnext.tw.rpi.edu/temporary/prefix/TemporaryGene")) ))
	real_gene_uris = retrieve_gene_dict(tmp_gene_uris)
	disease_file = open('disease.ttl', 'rw')
	data = disease_file.read()
	final_disease_file = open('disease_final.ttl', 'w')
	#Do an INSERT-DELETE for each temporary URI in the graph
	for tmp_gene_uri in tmp_gene_uris:

		print("TEMPORARY URI: " + tmp_gene_uri)
		real_gene_uri = ''
		if not real_gene_uris[tmp_gene_uri]:
			continue
		else:
			real_gene_uri = real_gene_uris[tmp_gene_uri].pop()
		print("REAL URI: " + real_gene_uri)
		symbol = tmp_gene_uri.split(TMP_PREFIX)[1]
		print("SYMBOL: " + symbol)
		data = data.replace("<%s>" % (tmp_gene_uri),
							 "<%s>" %(real_gene_uri))
		data = data.replace("TMP:%s" % (symbol),
							 "<%s>" %(real_gene_uri))
	final_disease_file.write(data)


"""
Will search for symbol in every filter listed in the 'filters' object.
Returns Entrez ID if available, HGNC ID as secondary when a match is found.
Otherwise returns None.
"""
def brute_force_biomart_for_entrez(symbol):
	# Taken from http://stackoverflow.com/questions/22270119/using-rpy2-and-biomart-in-django
	R.library("biomaRt")
	mart = R.useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

	filters = ['hgnc_symbol', 'ens_hs_gene', 'ens_hs_transcript', 'entrezgene', 'ensembl_gene_id', 
			   'clone_based_ensembl_gene_name', 'clone_based_ensembl_transcript_name', 
			   'clone_based_vega_gene_name', 'clone_based_vega_transcript_name',
			   'ensembl_gene_id']

	for bm_filter in filters:
		print("FILTER: " + bm_filter)
		tsv = R.getBM(attributes = StrVector(("hgnc_symbol", "hgnc_id","entrezgene")), 
						 filters = StrVector([bm_filter]), 
						 values = R.list(symbol), 
						 mart = mart)

		#Parse returned TSV and look to see if there's an Entrez ID;
		#if not then use the HGNC ID;
		#if that ain't workin' neither, then just use the bloody HGNC symbol

		print(tsv)
		try:
			if str(tsv[3][0]) != 'NA':
				return 'cortecon-neuralsci-org/cortecon/gene/' + str(tsv[3][0]) 
			else:
				pass
		except IndexError:
			pass
		try:
			if str(tsv[2][0]) != 'NA':
				return 'biomart-ensembl-org/biomart/hgnc/id/' + str(tsv[2][0]).split(':')[1]
			else: 
				pass
		except IndexError:
			pass
		try:
			if str(tsv[1][0]) != 'NA':
				return 'biomart-ensembl-org/biomart/hgnc/symbol/' + str(tsv[1][0]).split(':')[1] 
			else:
				pass
		except IndexError:
			pass
		try:
			if str(tsv[0][0]) != 'NA':
				return 'biomart-ensembl-org/biomart/ensembl/id/' + str(tsv[0][0]) 
			else:
				pass
		except IndexError:
			pass
	return None


'''
Create a dictionary of mappings from the provided symbol to a URI pulled from one of a variety of possible sources.
This will attempt to resolve the symbol against one of a number of sources in order to provide comprehensive coverage
of the genes present in the Setlr output.
'''
def retrieve_gene_dict(tmp_gene_uris):

	uris = {}
	for tmp_uri in tmp_gene_uris:
		#use regexp to extract symbol; reference against UMLS
		suffix = tmp_uri.split(TMP_PREFIX)[1]
		print("SYMBOL: " + suffix)

		#Check if the set already exists; if it does and a URI was found, just move on .
		#Otherwise instantiate an empty set for that temporary URI.
		if (tmp_uri in uris.keys()):
			if not(len(uris[tmp_uri]) == 0):
				continue
		else:
			uris[tmp_uri] = set()

		#First attempt to brute force biomart to find the Entrez ID
		ending = brute_force_biomart_for_entrez(suffix)
		if ending is not None:
			print("BIOMART HIT: " + str(ending))
			uris[tmp_uri] = set(['https://semnext.tw.rpi.edu/id/source/' + str(ending)])
			continue
		else:
		#Check Bio2RDF for the NCBI identifier first
			print("No BioMart hit, moving onto Bio2RDF...")
			try:
				uris[tmp_uri] = set( [res['@id'] for res in bio2rdf.bio2rdf_uri_for_symbol(suffix)] )
				print("BIO2RDF HIT: "+ str(uris[tmp_uri]))
			except EndPointInternalError:
				pass

		#Attempt a Bio2RDF free-text search as a last-ditch effort
		if (len(uris[tmp_uri]) == 0):
			print("Last-ditch Bio2RDf search...")
			try:
				#Create text variants that the string can match against
				variants = [suffix.lower(), suffix.replace('-', '').lower()]
				search_res = bio2rdf.search(suffix)
				bio2rdf_uri_to_semnext_uri = {
					'http://bio2rdf.org/ncbigene' : 'https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/gene/',
					'http://bio2rdf.org/hgnc' : 'https://semnext.tw.rpi.edu/id/source/bio2rdf-org/bio2rdf/hgnc/id/',
					'http://bio2rdf.org/hgnc.symbol': 'https://semnext.tw.rpi.edu/id/source/bio2rdf-org/bio2rdf/hgnc/symbol/'
				}

				#Iterate through the results; if a match against a valid variant of the gene symbol is found,
				#strip the value of interest from the URI and assign.
				for res in search_res:
					parts = res['@id'].split(':')
					print(res['@id'])
					print(parts)
					uri_prefix = bio2rdf_uri_to_semnext_uri[parts[0] + ":" + parts[1]]
					if res['value'].lower() in variants:
						print("BIO2RDF SEARCH HIT: " + uri_prefix + parts[2])
						uris[tmp_uri] = uri_prefix + parts[2] 
						continue
					else:
						pass
			except EndPointInternalError:
				pass
	
		#If nothing else has worked, then onto UMLS and lastly ReDrugS
		# if (len(uris[tmp_uri]) == 0):
		# 	uris[tmp_uri] = umls_db.lookup_gene_by_symbol(suffix)
		# if (len(uris[tmp_uri]) == 0):
		# 	uris[tmp_uri] = umls_db.lookup_gene_by_entrez_id(suffix)
		# if (len(uris[tmp_uri]) == 0):
		# 	uris[tmp_uri] = umls_db.lookup_gene_by_ensembl_id(suffix)
		# if (len(uris[tmp_uri]) == 0):
		# 	uris[tmp_uri] = redrugs.search(suffix)
		# 	print("REDRUGS URIs: " + str(uris[tmp_uri]))
		print(uris[tmp_uri])
	return uris


#TODO: remember to remove hard code from Setlr!!!
def main(argv=None):
	#Check arguments
	if argv is None:
		argv = sys.argv
		if argv is None:
			print("no arguments!")
			return 1

	#Call setlr on disease data
	# with open(argv[1]) as setl_file:
	# 	setl_graph = ConjunctiveGraph()
	# 	content = setl_file.read()
	# 	setl_graph.parse(data=content, format="turtle")
	# 	graphs = setlr._setl(setl_graph)

	#Check if Setl was successful
	'''
	if(os.path.isfile('./disease.ttl')):
		print("disease.ttl created") 
	else:
		print("disease.ttl not found")
		return 1
	'''
	#Replace temporary URIs
	g = Graph()
	g.parse(data=open('./disease.ttl').read(), format='turtle')
	g = replace_tmp_associated_genes(g)
	#g.serialize('./disease.ttl', format='trig')
	#for disease in g.triples(disease, None, None):

if __name__ == "__main__":
	sys.exit(main())
