# semnext-conversion
Raw conversion code for SemNExT knowledge graph

Description of Pipeline

1) Iterates through Jensen dataset and pulls in all genes listed for a DOID identifier (e.g. skips all ICD10 associations).
If the Entrez ID isn't available, uses the gene symbol instead (most commonly an HGNC symbol). This is our 'cleaned' Jensen data.

Note: at this point, there shouldn't be any genes that don't have an Entrez ID that we're interested in based on Nathan's data. 

2) Iterates through the Cortecon Gene Clock dataset and pulls all Entrez IDs, HGNC symbols from there.

3) Pulls dictionaries mapping DOID ID to disease name, 'see_also' URI.

4) Writes out every row of the cleaned Jensen data that has a matching Entrez ID and HGNC symbol from the Cortecon data;
includes the DOID ID, DOID name, evidence, source, confidence and see_also URI.

5) Invokes Setlr to transform this CSV into quads. Creates the following:

- quads for mapping https to http URIs (https://semnext.tw.rpi.edu/mapping/)
- subgraph for disease data (https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/disease/)
- subgraph for gene data (https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/gene/)

