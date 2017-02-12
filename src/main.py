#!/usr/bin/python
# -*- coding: utf-8 -*- 

"""
Generates an updated knowledge graph for data from the Cortecon and Jensen gene datasets
by expanding the existing N-Quad dump and replacing outdated assertions.
Accepts a target RDF file, two target CSV files to calculate the gene set intersection from.
Outputs an N-Quad file that captures the original information of the input RDF file,
plus new genes not formerly represented, minus outdated information.
"""

#Built-ins
import csv
import os
import sys
import logging
from shutil import copyfile
from pathlib import Path

#Local modules
import setlr

#RDFLib
from rdflib import ConjunctiveGraph, URIRef, Dataset, Graph
from rdflib.namespace import RDF

#rpy2
from rpy2.robjects import r as R
from rpy2.robjects.vectors import StrVector

#Set resource limits; it takes A LOT
#import resource
#resource.setrlimit(resource.RLIMIT_DATA, (2097152, 2621440)) #measured in kilobytes

#logging.basicConfig()


"""
Accepts a csv.DictReader; reads in Entrez IDs which are marked 'NULL'
and attempts to disambiguate them using the associated HGNC symbol.
Writes out Disease.database.table.cleaned.csv, which is the same file but with
the newly-found Entrez IDs.
Returns the number of rows that couldn't be cleaned.
"""
def clean_jensen_data(jensen_data):
	cleaned_output_headers = ['Entrez_ID','Gene_Symbol','Disease_Ontology_ID',
								'Confidence_Score','Source','Evidence', 'Type']
	cleaned_output_file = open('data/output/Disease.database.table.cleaned.csv', 'w')
	csvout = csv.DictWriter(cleaned_output_file, fieldnames=cleaned_output_headers)
	csvout.writeheader()

	#Iterate over jensen data. 
	#if the entrez ID is NULL, try to disambiguate.
	#If not a DOID or if Entrez ID not found, continue
	#TODO: uses Gene_Symbol for URI if Entrez ID not present; need to try to disambiguate in future
	count = 0
	for row in jensen_data:
		entrez = row['Entrez_ID']
		DOID = row['Disease_Ontology_ID']
		if entrez == 'NULL':
			#entrez = search_biomart_for_entrez(row['Gene_Symbol'])
			#if entrez is None:
			#	count += 1
			#	continue
			entrez = row['Gene_Symbol']
		if not 'DOID' in DOID:
			count += 1
			continue
		csvout.writerow({
			'Entrez_ID': entrez,
			'Gene_Symbol': row['Gene_Symbol'],
			'Disease_Ontology_ID': DOID,
			'Confidence_Score': row['Confidence_Score'],
			'Source': row['Source'],
			'Evidence': row['Evidence'],
			'Type': row['Type']
		})
	return count



"""
Will search for symbol in every filter listed in the 'filters' object.
Returns Entrez ID when a match is found.
Otherwise returns None.
NOT CURRENTLY USED, can be invoked later to try to disambiguate genes with missing Entrez IDs.
"""
def search_biomart_for_entrez(symbol):
	# Taken from http://stackoverflow.com/questions/22270119/using-rpy2-and-biomart-in-django
	R.library("biomaRt")
	mart = R.useMart(biomart="ensembl",dataset="hsapiens_gene_ensembl")

	filters = ['hgnc_symbol', 'ens_hs_gene', 'ens_hs_transcript', 'entrezgene', 'ensembl_gene_id', 
			   'clone_based_ensembl_gene_name', 'clone_based_ensembl_transcript_name', 
			   'clone_based_vega_gene_name', 'clone_based_vega_transcript_name',
			   'ensembl_gene_id']

	for bm_filter in filters:
		print("FILTER: " + bm_filter)
		tsv = R.getBM(attributes = StrVector(("entrezgene",)), 
						 filters = StrVector([bm_filter]), 
						 values = R.list(symbol), 
						 mart = mart)
		print(tsv)
		try:
			if str(tsv[0][0]) != 'NA':
				return str(tsv[0][0]) 
		except IndexError:
			pass
	return None



"""
Accepts CSV reader and returns a set of two-tuples with the Entrez ID and HGNC symbol for the gene.
"""
def get_cortecon_gene_superset(cortecon):
	genes = set()
	for row in cortecon:
		genes.add((row['Entrez_IDs'], row['Gene_Symbol']))
	return genes



"""
Writes 'cortecon-jensen-intersection.csv', which contains all genes present in both Cortecon,
plus associations in Jensen for genes contained in the Cortecon set.
"""
def write_intersection_dataset(cortecon_genes, cleaned_jensen_data, DOID_to_name, DOID_to_see_also):
	#Set up intersection file to be consumed by Setlr
	intersection_headers = ['Entrez_ID','Gene_Symbol','Disease_Ontology_ID',
							'Disease_Name', 'see_also',
							'Confidence_Score','Source','Evidence']
	intersection_file = open('data/output/cortecon-jensen-intersection.csv', 'w')
	intersection = csv.DictWriter(intersection_file, fieldnames=intersection_headers)
	intersection.writeheader()

	#Pull any genes that appear in both Cortecon and Jensen;
	#write intersection set out to new CSV file
	for row in cleaned_jensen_data:
		if (row['Entrez_ID'], row['Gene_Symbol']) in cortecon_genes:

			#Attempt to assign name and see_also since the DOID ID may not appear in disease_filtered.csv
			name = ''
			see_also = ''
			try:
				name = DOID_to_name[row['Disease_Ontology_ID']]
			except KeyError:
				pass
			try:
				see_also = DOID_to_see_also[row['Disease_Ontology_ID']]
			except KeyError:
				pass
			intersection.writerow( 
			{
				'Entrez_ID': row['Entrez_ID'], 
				'Gene_Symbol': row['Gene_Symbol'],
				'Disease_Ontology_ID': row['Disease_Ontology_ID'], 
				'Disease_Name': name,
				'see_also': see_also,
				'Confidence_Score': row['Confidence_Score'],
				'Source': row['Source'],
				'Evidence': row['Evidence']
			})
	intersection_file.close()



"""
Returns dictionary mapping tuples of the DOID symbol to the disease name
"""
def get_DOID_to_name(disease_filtered_file_path):
	DOID_to_name = {}
	with open(disease_filtered_file_path) as disease_filtered_file:
		disease_filtered_data = csv.DictReader(disease_filtered_file)
		for row in disease_filtered_data:
			DOID_to_name[row['DISEASE_ID']] = row['NAME']
	return DOID_to_name


"""
Returns dictionary mapping tuples of the DOID symbol to the see_also URI
"""
def get_DOID_to_see_also(disease_filtered_file_path):
	DOID_to_see_also = {}
	with open(disease_filtered_file_path) as disease_filtered_file:
		disease_filtered_data = csv.DictReader(disease_filtered_file)
		for row in disease_filtered_data:
			DOID_to_see_also[row['DISEASE_ID']] = row['ADDRESS']
	return DOID_to_see_also


"""
Invoke Setlr on the target .setl.ttl file. 
"""
def setl_data(setl_file_path):
	with open(setl_file_path) as setl_file:
			setl_graph = ConjunctiveGraph()
			content = setl_file.read()
			setl_graph.parse(data=content, format='turtle')
			graphs = setlr._setl(setl_graph)
	return 0



"""
Remove all quads in the cortecon/disease and cortecon/gene subgraphs of the original 
semnext_dump.nq; insert genes into cortecon/gene and diseases into cortecon/disease. 
"""
def update_knowledge_graph():
	print("Beginning graph update.")

	#Get quads we don't reproduce but want to keep
	preserve_desired_cortecon_graph_quads()
	print("Preserved quads extracted.")

	#Clean out anything that is in a cortecon/disease or cortecon/gene graph
	wipe_cortecon_graphs()
	print("Knowledge graph cleaned.")

	#merge cleaned and preserved quads; write them out!
	merge_cleaned_and_preserved_quads()
	print("Cleaned and preserved quads merged.")

	#merge the cleaned/preserved quads with the new data
	merge_cleaned_preserved_quads_and_disease_data()
	print("Final graph merged; success!")

"""
Accepts an RDF Graph of quads; returns any quads coming from the
https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon graph that
matches against the predicate list. This is done in order to preserve data from the cluster
analysis that is not generated as part of this update.
It's also currently unclear what TCONS genes are, so I'm going to preserve all of the information for the
meantime in the event that it's actually rather important.
"""
#TODO: this is disgustingly inefficient, should trim this down in the future.
def preserve_desired_cortecon_graph_quads():

	#Check if file already present
	preserved_quads_path = Path("data/output/preserved-quads.nq")
	if preserved_quads_path.is_file():
		return

	#TODO: Read in data as stream, don't load in like this because it's wicked inefficient
	#Load up SemNExT dump file
	classes = [
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#FuzzyClusterMembership'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#PrincipalComponentAnalysis'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#PCAOutput'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#PCAVector'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#PCAScore'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#Cluster'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#RelativeActivation'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#Experiment')
	]

	predicates = [
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#fuzzyMemberOf'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#coordinate'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#hasAnalysis'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#inCluster'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#score'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#frequency'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#durationSinceStart'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#z-score'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#dimension'),
		URIRef('http://semnext.tw.rpi.edu/ontology/semnext#memberOf'),
		URIRef('http://open.vocab.org/terms/subjectDiscriminator'),
		URIRef('http://purl.org/dc/terms/isReferencedBy'),
		URIRef('http://rdfs.org/ns/void#inDataset'),
		URIRef('http://purl.org/dc/terms/identifier')
	]

	#Add in anything that falls under the SemNExT class structure for cluster analyses;
	#add any properties attached to the instance
	preserved_quads = Graph()
	with open('data/quads/semnext-dump.nq') as fin:
		fout = open('data/output/preserved-quads.nq','w')
		for line in fin:
			quad = Dataset()
			quad.parse(data=line, format='nquads')

			#Match against our instance type declarations
			for s, p, o, context in quad.quads((None, RDF.type, None,
				'https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/')):
				if o in classes:
					fout.write(line)

			#look for the predicates we're interested in, plus capture anything with TCON in the name
			for s, p, o, context in quad.quads((None, None, None,
				'https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/')):
				if p in predicates or 'TCONS' in s or 'TCONS' in o:
					fout.write(line)
		fout.close()


"""
Once the desired quads have been preserved, eliminate all quads in the knowledge graph from
the cortecon disease and gene subgraphs.
"""
def wipe_cortecon_graphs():
	#Check if file already present
	cleaned_quads_path = Path("data/output/cleaned-quads.nq")
	if cleaned_quads_path.is_file():
		return

	# knowledge_graph = Dataset()
	# knowledge_graph.parse(format='nquads', source='data/quads/semnext-dump.nq')
	# wiped_knowledge_graph = Dataset()
	prefix = 'https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/'
	with open('data/quads/semnext-dump.nq') as fin:
		fout = open('data/output/cleaned-quads.nq', 'w')
		for line in fin:
			quad = Dataset()
			quad.parse(data=line, format='nquads')
			for s, p, o, context in quad.quads((None, None, None, None)):
				if (prefix in context.n3()):
					continue
				else:
					fout.write(line)
		fout.close()



def merge_cleaned_and_preserved_quads():
	#Check if file already present
	merged_preserved_cleaned_quads_path = Path("data/output/merged-preserved-cleaned-quads.nq")
	if merged_preserved_cleaned_quads_path.is_file():
		return

	#Copy the cleaned quads file
	copyfile('data/output/cleaned-quads.nq', 'data/output/merged-preserved-cleaned-quads.nq')
	with open('data/output/preserved-quads.nq') as fin:
		with open('data/output/merged-preserved-cleaned-quads.nq', 'a') as fout:
			for line in fin:
				fout.write(line)


def merge_cleaned_preserved_quads_and_disease_data():
	#make a copy of the merged-preserved quads file
	copyfile('data/output/merged-preserved-cleaned-quads.nq', 'data/output/semnext-graph.nq')

	disease_trig = Dataset()
	disease_trig.parse(format='trig', source='data/output/disease.trig')
	print("disease.trig loaded.")

	#Iterate through the merged and preserved quads; if the quad appears in the disease data,
	#drop it from the disease data so that 
	with open('data/output/semnext-graph.nq', 'r') as final_graph:
		for line in final_graph:
			quad = Dataset()
			quad.parse(data=line, format='nquads')
			for s, p, o, context in quad.quads((None, None, None, None)):
				if ((s, p, o, context) in disease_trig):
					disease_trig.remove((s, p, o, context))
	print("First pass success; attempting merge...")

	#Finally, dump all of the N-Quads from the cleaned disease data into the merged dataset
	with open('data/output/semnext-graph.nq', 'a') as final_graph:
		data = disease_trig.serialize(format='nquads')
		print("disease_trig serialized, appending...")
		final_graph.write(data)

	#We did it!


"""
Main method expects:
	argv[1]: path to .setl.ttl file that describes ingest of the intersection dataset
	argv[2](optional): path to RDF file to modify. If provided, will strip any existing information
						about each entity provided in the Setlr output and replace it with the contents of the output.
"""
def main(argv=sys.argv):
	#Pull in datasets
	gene_clock_data_headers = ['Entrez_IDs', 'Gene_Symbol' ,
								'd0','d7','d12','d19','d26','d33','d49','d63','d77',
								'Cluster','Radius','Angle',
								'Autism','Holoprosencephaly','Microcephaly','Lissencephaly','Alzheimers','Tauopathy']
	jensen_data_headers = 	  ['','Entrez_ID','Gene_Symbol','Disease_Ontology_ID',
								'Confidence_Score','Source','Evidence', 'Type']						
	gene_clock_data = csv.DictReader(open('data/tabular/GeneClockData.csv'))
	jensen_data = csv.DictReader(open('data/tabular/Restricted.disease.table.csv'))

	#Clean jensen data of NULL values for Entrez IDs; remove any rows without DOID ID;
	#Write out to 'Disease.database.table.cleaned.csv'
	miss_count = clean_jensen_data(jensen_data)
	print("Jensen data cleaned; " + str(miss_count) + " rows could not be cleaned.")

	#Retrieve set of Entrez IDs appearing in Cortecon
	cortecon_genes = get_cortecon_gene_superset(gene_clock_data)

	#Retrieve dictionaries mapping DOID to names, to the associated jensendata.com URI for rdf:see_also
	disease_filtered_file_path = 'data/tabular/disease_filtered.csv'
	DOID_to_see_also = get_DOID_to_see_also(disease_filtered_file_path)
	DOID_to_name = get_DOID_to_name(disease_filtered_file_path)

	#Write out intersection of Cortecon and Jensen datasets, plus disease name & see_also
	#From disease_filtered.csv 
	cleaned_jensen_data = csv.DictReader(open('data/output/Disease.database.table.cleaned.csv'))
	print("Jensen data read.")
	write_intersection_dataset(cortecon_genes, cleaned_jensen_data, DOID_to_name, DOID_to_see_also)

	#Invoke Setlr on intersection file
	#setl_data(argv[1])

	#Replace contents of RDF file if provided and write out updated graph
	update_knowledge_graph()

	return 0

if __name__ == '__main__':
	sys.exit(main())
