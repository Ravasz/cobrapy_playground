'''
Created on 1 Feb 2019

@author: mate

- take in uniprot IDs from processed proteinGroups file.
- convert mouse uniprot IDs to mouse gene symbols
- convert mouse gene symbols to human gene symbols
- convert mouse gene symbols with no human counterpart by blast-ing the sequence for a human homologue
- convert human gene names to HGNC IDs
'''
