cosmic2vcf
--------------
Convert COSMIC structural variant TSV files to VCF format  
Eric T Dawson   
June 2019

### Intro
The Catalog of Somatic Mutations in Cancer contains a publicly-available
set of variant calls seen across many tumors. Some of these files (specifically
coding and non-coding small variants) are in the standard VCF format while others,
most notably the structural variant calls, are in a non-standard TSV format. This
frustrates their usage as inputs to other programs.

This repo contains a python script to convert the CosmicStructExport.tsv file
to a standard VCF4.3 file that has the required SV info field tags such as 
SVLEN, SVTYPE, END, SPAN, etc. This VCF file can then be sorted/bgzipped/indexed
and used in downstream analyses.

### Usage
```
python cosmic_structexport_to_vcf.py -i CosmicStructExport.tsv > csv.vcf
```

Currently, this script handles roughly 3/4 of the variants in COSMIC. It does no
special handling of fold-back inversions or any other variants.

### License
MIT
