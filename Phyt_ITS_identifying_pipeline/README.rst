READme info Phyt_ITS_identifying_pipeline
======================================================
author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

basic usage:

run an individual shell script, which acts as a config file, that then call the main 
script below 

``./Phy_CONFIG_file.sh``

``./Phyt_ITS_identify_pipline.sh``

Important parameters to be set:
some variables in here need to be set by the user.
see top section marked: ## THESE VARIABLE NEED TO BE FILLED IN BY USER !!!
until: # NOTHING TO BE FILLED IN BY USER FROM HERE!



why?: We are interested in identifying Phyt species based on the ITS
regions. Phy species could become more previlent in the UK due to climate change and
the importation of plants.

How?: Following plant sampling, DNA extraction, primers amplify the ITS region.
PCR amplified products are Illumina sequenced. The reads are QC, trmimmed, QC,
assembled together (reads are overlapping) and clustered with a "known" database of Phyt
ITS regions. The resulting clusters are used to identify the Phy species present.

Results?: NOT FINISHED. Basically, the clustering is one cluster per line, tab separated.



``Requires:``
	1) Internet access for $ wget
	2) Trimmomatic (read quality trimming) http://www.usadellab.org/cms/?page=trimmomatic
	3) Fastqc http://www.bioinformatics.babraham.ac.uk/projects/download.html
	4) python2.7 or greater. Python 3.5 is recommended. https://www.python.org/downloads/
	5) PEAR (assemble overlapping reads) https://github.com/xflouris/PEAR  http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz 
	6) Swarm (clustering) https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf
	7) To draw graphs: No essential but Numpy, Matplotlib, Scipy
	8) BLESS error correction (difficult to install) - Not essential
	9) SPAdes alternative error correction. Much easier to install .
		



