#!/bin/bash
# put out files in current directory
#$ -cwd
#$ -l hostname="n13*"
# Deliver notifications to the following address
# Send notifications when the job begins (b), ends (e), or is aborted (a)
#dollar -M email.address@somewhere .... currently doesnt work
#dollar -m a 
#Abort on any error,
#set -euo pipefail
#(and thus quit the script right away)

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile

cd ~/scratch/tree_health/ITS_ratio
export TMP=~/scratch/${USER}_${JOB_ID}

##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

species=Phytophthora_cambivora

genome_prefix=GCA_000443045.1_PhytoCambi_1.0

genome_fasta=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000443045.1_PhytoCambi_1.0/GCA_000443045.1_PhytoCambi_1.0_genomic.fna.gz

genome_GFF=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000443045.1_PhytoCambi_1.0/GCA_000443045.1_PhytoCambi_1.0_genomic.gff.gz

read_1_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR870/SRR870485/SRR870485_1.fastq.gz

read_2_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR870/SRR870485/SRR870485_2.fastq.gz

trimmomatic_path=~/Downloads/Trimmomatic-0.32

SRA_prefix=SRR870485

repository_path=~/misc_python/THAPBI/ITS_region_genomic_coverage

num_threads=12


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################

export species
export genome_prefix
export genome_fasta
export genome_GFF
export read_1_link
export read_2_link
export trimmomatic_path
export SRA_prefix
export repository_path
export num_threads


./ITS_genomic_coverage_ratio_finding.sh
