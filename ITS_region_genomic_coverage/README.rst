READme info ITS_region_genomic_coverage
======================================================
author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

basic usage:

``./P.nicotiana_ITS_ratio_finding.sh`` 

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
Internet access for $ wget
Trimmomatic (read quality trimming)
BLAST
python
bowtie2
samtools
bedtools 
(for poorly formamted GFF files - try reformmatting with GenomeTools)
