# semnext-conversion
Raw conversion code for SemNExT knowledge graph

<h2>First order of business</h2>

You'll need the current version of the SemNExT knowledge graph in order to perform the update. This isn't included in the repo since the graph is exceptionally large.

Run the following commands:

```
$ cd ./data/quads
$ curl https://semnext.tw.rpi.edu/data/source/semnext-tw-rpi-edu/semnext/version/${VERSION}/semnext-${VERSION}.nq.gz" > semnext-dump.nq.gz
$ gunzip ./semnext-dump.nq.gz
```
The `${VERSION}` variable should be replaced with the latest version of the knowledge graph. As of 02/12/2017, this is `0.1.3`. 

<h2>Description of Pipeline</h2>

1) Iterates through Jensen dataset and pulls in all genes listed for a DOID identifier (e.g. skips all ICD10 associations).
If the Entrez ID isn't available, uses the gene symbol instead (most commonly an HGNC symbol). This is our 'cleaned' Jensen data.

- `TODO`: at this point, there shouldn't be any genes that don't have an Entrez ID that we're interested in based on Nathan's data. Probably a good idea to modify the script so that we don't accidentally use a gene symbol anywhere in a URI. 

2) Iterates through the Cortecon Gene Clock dataset and pulls all Entrez IDs, HGNC symbols from there.

3) Pulls dictionaries mapping DOID ID to disease name, 'see_also' URI.

4) Writes out every row of the cleaned Jensen data that has a matching Entrez ID and HGNC symbol from the Cortecon data;
includes the DOID ID, DOID name, evidence, source, confidence and see_also URI.

5) Invokes Setlr to transform this CSV into quads. Creates the following in `disease.trig`:

- quads for mapping https to http URIs (https://semnext.tw.rpi.edu/mapping/)
- subgraph for disease data (https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/disease/)
- subgraph for gene data (https://semnext.tw.rpi.edu/id/source/cortecon-neuralsci-org/cortecon/gene/)
- `TODO`: This step needs to be done manually. The `main.py` script bombs out halfway through, at which point you invoke Setlr, comment out the actual call to Setlr in `main.py` and re-run the whole thing. This isn't terrible inefficient (adds ~30 seconds in total) but it is confusing and doesn't lend itself to integration. MUST FIX.

6) Begins update of knowledge graph. This is where it gets hairy in the code since everything needs to be streamed to avoid hitting memory limits (`TODO`: memory consumption still goes up during the final merge, not sure why. See if we can fix this; would reduce needed memory by ~1Gb; otherwise you need to bump the settings in the Vagrantfile).

- `Note`: if you update the tabular files, you need to delete any existing files in `output` before you run again; otherwise, the script will see them and skip them (this was done to reduce run times when iterating on existing work).

- Extract any quads from Cortecon subgraphs in the existing SemNExT dump that pertain to the cluster analyses or TCONS data, write these out to `preserved-quads.nq`
- Strip the SemNExT dump of all Cortecon subgraphs (there are currently subgraphs for each and every disease and gene); write this stripped dump to `cleaned-quads.nq`
- Merge `preserved-quads.nq` and `cleaned-quads.nq` and write to `merged-preserved-cleaned-quads.nq`
- Merge `merged-preserved-cleaned-quads.nq` and `disease.trig`, write to the final `semnext-graph.nq`

Voila, a complete knowledge graph. Now gzip this and drop it at `SemNExT/data/source/semnext-graph.nq.gz` and start up your Docker service.
