# Norovirus
Nextstrain analysis of Norovirus

## Introduction
This is the Nextstrain build for Norovirus

## Data Curation
All sequence data is from Vipr or Genbank. The full Norovirus genomic length is ~7,547 bp long. In this build, we filtered for human Norovirus sequences that are at least 5032bp long (2/3 of the full length). We ended up with a dataset of 2734 sequences from 1968-2022, from 42 countries.
### Obtaining large dataset metadata and sequence output files
1. parse data/GenomicFastaResults.fasta using augur parse --sequences data/GenomicFastaResults.fasta --output-sequences results/sequence_output.fasta --output-metadata results/metadata.tsv --fix-dates monthfirst --fields {"strain","strain_name","segment","date","host","country"} 
2. convert resulting metadata.tsv file to NDJSON format.
3. fix the dates using the transform-dates-field script by running ./transform-date-fields.py --expected-date-formats "%Y_%m_%d" "%Y_%m_%dT%H:%M:%SZ" "%Y_%m" --date-fields date < results/metadata.NDJSON.
4. convert resulting metadata.NDJSON file to tsv file, naming it metadata_parsed_big.tsv
5. the sequence_output.fasta generated in this step is unnecessary. We will be using the sequence_output.fasta from the next step, so we can delete this file.
### Obtaining smaller dataset sequence output file 
1. parse data/sequence.fasta using augur parse --sequences data/GenomicFastaResults.fasta --output-sequences results/sequence_output.fasta --output-metadata results/metadata.tsv --fix-dates monthfirst --fields {"strain","strain_name","segment","date","host","country"} 
2. the resulting metadata file from the parse command is unnecessary, we will be obtaining that from the third step
### Joining the datasets
1. run the join.py script using python join.py. There were redundant columns as a result, so I deleted the redundant ones in excel.
2. the sequence_output file had some sequences that did not show up in the metadata. I had to go in and manually remove them. There must be a better way to do this. 

## Further Reading
Relevant papers for further reading:
* [Norwalk Virus Minor Capsid Protein VP2 Associates within the VP1 Shell Domain](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624303/)
* [Deep Sequencing of Norovirus Genomes Defines Evolutionary Patterns in an Urban Tropical Setting](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4178781/)
