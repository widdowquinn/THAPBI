READme info Phyt_ITS_identifying_pipeline
======================================================
author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

basic usage:

run an individual shell script, which acts as a config file, that then call the main 
script below 
# config file - see bottom of here for how to generate this file. 
``./Phy_CONFIG_file.sh``

#which calls
``./Phyt_ITS_identify_pipline.sh``

or

``./Phyt_ITS_identify_pipline_with_error_correction.sh``

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
	1) Trimmomatic (read quality trimming) http://www.usadellab.org/cms/?page=trimmomatic
	2) Fastqc http://www.bioinformatics.babraham.ac.uk/projects/download.html
	3) python2.7 or greater. Python 3.5 is recommended. https://www.python.org/downloads/
	4) PEAR (assemble overlapping reads) https://github.com/xflouris/PEAR  http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz 
	5) Swarm (clustering) https://github.com/torognes/swarm/blob/master/man/swarm_manual.pdf
	6) To draw graphs: No essential but Numpy, Matplotlib, Scipy... pip install module_name
	7) BLESS error correction (difficult to install) - Not essential. NOT CURRENTLY USED
	8) SPAdes alternative error correction. Much easier to install . http://bioinf.spbau.ru/en/content/spades-download-0 
		

generate config file
=====================
Either use the example config file as a template and fill in your experiment details.
Dont forget to make the file executable:
$chmod a+x names_of_config_file.sh``


or

use the python script:

``python generating_config_files.py``



