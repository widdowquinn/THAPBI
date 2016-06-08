READme info ITS_region_genomic_coverage
======================================================
author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

basic usage:

``./ITS_genomic_coverage_ratio_finding.sh`` 

All regions in this file between ####
need to be filled in by the user.

why?: We are interested in the raio of coverage of all identified
ITS regions in a genome versus "normal" genes - whatever they are.
Normal genes, as these cant be defined, will be treated as all other
genes, except RNA annotated genes (if annotated in GFF).

How?: Illumina reads are mapped back to the genome. Using the GFF3
gene coordingates, the coverage in obtained. A GFF3 file is generated
based on the BLAST results (ITS versus genome). The coordinates 
are used to obtain the genomic read coverage for these ITS regions.

Results?: The mean number of reads the map to "normal genes" versus 
ITS can be compared. This ratio give us an indication of the total
number of ITS regions that should be there, assuming genomic read 
coverage is normal and representative of gene number.


``Requires:``
	1) Internet access for $ wget
	2) Trimmomatic (read quality trimming)
	3) BLAST
	4) python
	5) bowtie2
	6) samtools
	7) bedtools 
	7b)(for poorly formamted GFF files - try reformmatting with GenomeTools)
