#!/usr/bin/env bash
#
# Title: pipline to identify Phy species by clustering with ITS regions
# Author:  Peter Thorpe
# (C) The James Hutton Institute 2016


##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Name_of_project=testing

read_name_prefix="M01157:20:000000000-D07KA:1:"

left_read_file=DNAMIX_S95_L001_R1_001.fastq.gz

right_read_file=DNAMIX_S95_L001_R2_001.fastq.gz

trimmomatic_path=~/Downloads/Trimmomatic-0.32

repository_path=~/misc_python/THAPBI/Phyt_ITS_identifying_pipeline

num_threads=4


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################
#not needed but, may use it as a configuration style script later
export trimmomatic_path
export repository_path
export num_threads

cd $HOME/misc_python/THAPBI/Phyt_ITS_identifying_pipeline/data

# Create directory for output
mkdir fastq-join_joined

# QC "forward" and "reverse" reads
echo "running fastqc on raw reads"
cmd_fastqc="fastqc ${left_read_file}" 
echo ${cmd_fastqc}
eval ${cmd_fastqc}
cmd_fastqc2="fastqc ${right_read_file}" 
echo ${cmd_fastqc2}
eval ${cmd_fastqc2}
echo "cmd_fastqc done"

#quality trim the reads
echo "Trimming:"
echo "have you added to the adapter seq to the triming database file?"
echo "do we want to keep the adapters in right until the end to track the seq origin?"
cmd_trimming="java -jar ${trimmomatic_path}/trimmomatic-0.32.jar PE -threads ${num_threads} -phred33 ${left_read_file} ${right_read_file} ${Name_of_project}_R1.fq.gz unpaired_R1.fq.gz ${Name_of_project}_R2.fq.gz unpaired_R2.fq.gz ILLUMINACLIP:${repository_path}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75" 
echo ${cmd_trimming}
eval ${cmd_trimming}
echo "Trimming done"

#remove unpaired reads
echo "removing unpaired files."
rm unpaired*

# Conduct QC on trimmed reads
echo "running fastqc on trimmed reads"
cmd_fastqc="fastqc ${Name_of_project}_R1.fq.gz" 
echo ${cmd_fastqc}
eval ${cmd_fastqc}
cmd_fastqc2="fastqc ${Name_of_project}_R2.fq.gz" 
echo ${cmd_fastqc2}
eval ${cmd_fastqc2}
echo "cmd_fastqc done"

# Join trimmed paired-end read files together
# use a program called PEAR
echo "puttin the reads together:"
cmd_pear="pear -f ${Name_of_project}_R1.fq.gz -r ${Name_of_project}_R2.fq.gz -o ${Name_of_project}_PEAR" 
echo ${cmd_pear}
eval ${cmd_pear}
echo "cmd_pear done"

# Move the joined read files to the output subdirectory
mv ${Name_of_project}_PEAR.assembled.fastq ${Name_of_project}_outfiles

# Change directory to the output subdirectory
cd ${Name_of_project}_outfiles

# Conduct QC on joined reads
fastqc fastqjoin.join.fastq

# Convert read format from FASTQ to FASTA
# convert_format comes from seq_crumbs: https://github.com/JoseBlanca/seq_crumbs
echo "converting the fq to fa file"
cmd_convert="convert_format -t fastq -o fastqjoin.join.fasta -f fasta ${Name_of_project}_PEAR.assembled.fastq" 
echo ${cmd_convert}
eval ${cmd_convert}
echo "cmd_convert done"


#cluster with swarm.
#need to rename the reads.

echo "renaming the fa id. Swarm doesnt like these as they are"
cmd_rename="python ${repository_path}/re_format_fasta.py -f fastqjoin.join.fasta --prefix ${read_name_prefix}" 
echo ${cmd_rename}
eval ${cmd_rename}
echo "cmd_rename done"

echo "running swarm: swarm [OPTIONS] [filename]"
cmd_swarm="swarm -t ${num_threads} -d 1 -o swarm_results fastqjoin.join_alt.fasta" 
echo ${cmd_swarm}
eval ${cmd_swarm}
echo "cmd_swarm done"
 
swarm -t ${num_threads} -d 1 -o Phy_ITSregions_all_20160601.fixed02.fasta.swarm_clusters Phy_ITSregions_all_20160601.fixed02.fasta 


python ${repository_path}/parse_clusters.py -i ${repository_path}/data/fastq-join_joined/Phy_ITSregions_all_20160601.fixed02.fasta.swarm_clusters

	   
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
python ~/bin/pick_closed_reference_otus.py -i fastqjoin.join.fasta -r reference.fasta -o otus/
python ~/bin/pick_otus.py  -i fastqjoin.join.fasta -r reference.fasta -m uclust_ref -s 0.99
