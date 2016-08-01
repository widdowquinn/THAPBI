#!/usr/bin/env bash
#
# Title: pipline to identify Phy species by clustering with ITS regions
# Author:  Peter Thorpe
# (C) The James Hutton Institute 2016


##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

trimmomatic_path=~/Downloads/Trimmomatic-0.32

repository_path=~/misc_python/THAPBI/ITS_region_genomic_coverage

num_threads=4


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################

export trimmomatic_path
export repository_path
export num_threads



# Create directory for output
mkdir fastq-join_joined

# Copy forward and reverse reads into file with 'standard' naming

cp *R1* temp_r1.fastq.gz
cp *R2* temp_r2.fastq.gz

# Conduct QC on collated "forward" and "reverse" reads
fastqc temp_r1.fastq.gz
fastqc temp_r2.fastq.gz

#quality trim the reads
echo "Trimming:"
echo "have you added to the adapter seq to the triming database file?"
cmd_trimming="java -jar ${trimmomatic_path}/trimmomatic-0.32.jar PE -threads ${num_threads} -phred33 temp_r1.fastq.gz temp_r2.fastq.gz R1.fq.gz unpaired_R1.fq.gz R2.fq.gz unpaired_R2.fq.gz ILLUMINACLIP:${repository_path}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75" 
echo ${cmd_trimming}
eval ${cmd_trimming}
echo "Trimming done"

#remove unpaired reads
rm unpaired*

# Conduct QC on trimmed reads
fastqc R1.fq.gz
fastqc R2.fq.gz



# Move untrimmed, collated read files to output subdirectory
mv DCM1.fastq fastq-join_joined
mv DCM2.fastq fastq-join_joined

# Join trimmed paired-end read files together
# use a program called PEAR
echo "puttin the reads together:"
cmd_pear="pear -f R1.fq.gz -r R2.fq.gz -o fastq-join_joined" 
echo ${cmd_pear}
eval ${cmd_pear}
echo "cmd_pear done"

# Move the joined read files to the output subdirectory
mv fastq-join_joined.assembled.fastq fastq-join_joined

# Remove all .zip files in current directory
rm *.zip

# Change directory to the output subdirectory
cd fastq-join_joined

# Conduct QC on joined reads
fastqc fastqjoin.join.fastq

# Convert read format from FASTQ to FASTA
# convert_format comes from seq_crumbs: https://github.com/JoseBlanca/seq_crumbs
convert_format -t fastq -o fastqjoin.join.fasta -f fasta fastq-join_joined.assembled.fastq

# Trim something?
# This seems to be adaptor sequences, from the Materials and Methods
python trim_longitudes.py # this script take off the first and last 20 bp? needed???

# Cluster FASTA reads with BLASTCLUST and [convert output to FASTA]?
#blastclust -L 0.90 -S 99 -i fastqjoin.join.fasta -o fastqjoin.join.blastclust99.lst -a 4 -p F


#cluster with swarm.
#need to rename the reads.
python 
swarm [OPTIONS] [filename]

 swarm -d 1 -o swarm_results test.fasta
swarm -t ${num_threads} 
		   
python blastclust_lst2fasta.py

# Create subdirectories for output
mkdir data
mkdir scripts
mv *.py scripts
mkdir align

# Move output to subdirectories
mv *OTU* align
mv *fastq data
mv *lst data
mc *.fasta data
mc *txt data
rm *.zip
mv *error* data

# Change to align subdirectory, and run MUSCLE alignment
cd align
for i in fasta
do
  muscle -in $i -out $i"_muscle.fasta"
done

# Use QIIME to pick OTUs
pick_closed_reference_otus.py -i fastqjoin.join.fasta -r reference.fasta -o otus/
pick_otus.py -i fastqjoin.join.fasta -r reference.fasta -m uclust_ref -s 0.99
