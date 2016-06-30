#!/bin/bash
# put out files in current directory
#$ -cwd
#$ -l hostname="n13*"
# Deliver notifications to the following address
# Send notifications when the job begins (b), ends (e), or is aborted (a)
#dollar -M email.address@somewhere .... currently doesnt work
#dollar -m a 
#Abort on any error,
set -euo pipefail
#(and thus quit the script right away)

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile

cd ~/scratch/tree_health/ITS_ratio
export TMP=~/scratch/${USER}_${JOB_ID}

##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

species=Phytophthora_cinnamomi

genome_prefix=GCA_001314365.1_MP94-48v2_genomic

genome_fasta=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_001314365.1_MP94-48v2/GCA_001314365.1_MP94-48v2_genomic.fna.gz

genome_GFF=ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_001314365.1_MP94-48v2/GCA_001314365.1_MP94-48v2_genomic.gff.gz

read_1_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR097/SRR097634/SRR097634_1.fastq.gz

read_2_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR097/SRR097634/SRR097634_2.fastq.gz

trimmomatic_path=~/Downloads/Trimmomatic-0.32

SRA_prefix=SRR097634

path_to_ITS_clipping_file=~/misc_python/THAPBI/ITS_region_genomic_coverage

num_threads=12


# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################

#put these variable to file for logging purposes.
echo "USER VARIABLES ARE: 
genome_prefix = ${genome_prefix}, 
read_1_link = ${read_1_link}, 
read_2_link = ${read_2_link}, 
genome_fasta = ${genome_fasta} 
genome_GFF = ${genome_GFF}
num_threads = ${num_threads}
path_to_ITS_clipping_file = ${path_to_ITS_clipping_file}
trimmomatic_path = ${trimmomatic_path}
"

mkdir ${genome_prefix}
cd ${genome_prefix}


#reads genome and GFF, gunzip 

# i get my reads from: http://www.ebi.ac.uk/ena/data/warehouse/search
# use the advance.  Species ENTER, Lib.layout=paired,  platform=illumina, lib.source=genomic > search
echo STEP1: Downloading reads and genome files
wget ${genome_fasta}
wget ${genome_GFF}
gunzip *.gz

wget ${read_1_link}
wget ${read_2_link}

# EXAMPLE: Phytophthora_kernoviae.GCA_000333075.1.31.dna.genome.fa.gz => Phytophthora_kernoviae.GCA_000333075.1.31.
#gunzip ${genome_prefix}*


# blast to get representative ITS regions.

echo " STEP2: blast searches"
cmd="makeblastdb -in ${genome_prefix}*.fa -dbtype nucl" 
echo ${cmd}
eval ${cmd}

cmd2="blastn -query ${path_to_ITS_clipping_file}/Phy_ITSregions_all_20160601.fasta -db ${genome_prefix}*.fna -outfmt 6 -out n.Pi_ITS_vs_${genome_prefix}.out" 
echo ${cmd2}
eval ${cmd2}


# prepare ITS gff. python:
echo "STEP 3: prepare ITS gff. python
"
cmd_python_ITS="python ${path_to_ITS_clipping_file}/generate_ITS_GFF.py --blast n.Pi_ITS_vs_${genome_prefix}.out --prefix ${genome_prefix} -o ${genome_prefix}.ITS.GFF" 
echo ${cmd_python_ITS}
eval ${cmd_python_ITS}

cat ${genome_prefix}.ITS.GFF | uniq | sort -k1,1 -k4,4 -k5,5 > temp.out
mv temp.out ${genome_prefix}.ITS.GFF

# genearate consensus blast hit ITS GFF file from the sorted file above.
wait
cmd_python_ITS_consensus="python ${path_to_ITS_clipping_file}/filter_GFF.py --gff ${genome_prefix}.ITS.GFF -o ${genome_prefix}.ITS.consensus.GFF"
echo ${cmd_python_ITS_consensus}
eval ${cmd_python_ITS_consensus}


wait
#quality trim the reads
echo "Trimming:"
cmd_trimming="java -jar ${trimmomatic_path}/trimmomatic-0.32.jar PE -threads ${num_threads} -phred33 ${SRA_prefix}_1.fastq.gz ${SRA_prefix}_2.fastq.gz R1.fq.gz unpaired_R1.fq.gz R2.fq.gz unpaired_R2.fq.gz ILLUMINACLIP:${path_to_ITS_clipping_file}/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:45" 
echo ${cmd_trimming}
eval ${cmd_trimming}
echo "Trimming done"

wait
#remove the unpaired reads....
rm  unpaired_R*

#map to the genome
mkdir $TMP

# index genome
echo "index genome"
cmd_index="bowtie2-build --quiet -f ${genome_prefix}*.fna bowtie_index_files" 
echo ${cmd_index}
eval ${cmd_index}

# randomly assign mutliple mapping reads ...  http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-p
#OLD command: bowtie2 --very-sensitive --non-deterministic --seed 1 --no-mixed --no-unal -p 8 -x Pi -1 R1.fq.gz -2 R2.fq.gz -S P.nicotiana.sam

#pipe SAM output on stdout into samtools to make it into BAM
#TODO: Put the temp unsorted BAM file on the local hard drive scratch space under /mnt/scratch
# index genome
echo "MAPPING"
cmd_mapping="bowtie2 --very-sensitive --non-deterministic --seed 1 --no-mixed --no-unal -p ${num_threads} -x bowtie_index_files -1 R1.fq.gz -2 R2.fq.gz | samtools view -@ ${num_threads} -S -b -o $TMP/tmp_unsorted.bam -" 
echo ${cmd_mapping}
eval ${cmd_mapping}
echo "mapping fininshed"

# convert to sorted bam.
#samtools view -@ 8 -S -b -o tmp_unsorted.bam P.nicotiana.sam
echo "sort bam file"
cmd_sort="samtools sort -@ ${num_threads} $TMP/tmp_unsorted.bam ${genome_prefix}"
echo ${cmd_sort}
eval ${cmd_sort}

#index bam file
samtools index ${genome_prefix}.bam


# Clean up, doing this last in case fails
rm $TMP/tmp_unsorted.bam
# WARNING - Blindly deleting folders based on an encironemnt variable is DANGEROUS!
# There I am NOT doing: r m -r f $ T M P
# rmdir will only work if the directory is emtpy (and it should be)
rmdir $TMP

# get only the genes, not bothered about other stuff ...
echo "prepare GFF for genes only"
cmd_python_gene_to_gff=" python ${path_to_ITS_clipping_file}/THAPBI/ITS_region_genomic_coverage/get_genes_from_GFF.py --gff ${genome_prefix}*gff3 -o ${genome_prefix}.gene.gff"
echo ${cmd_python_gene_to_gff}
eval ${cmd_python_gene_to_gff}
 
#old commands - doesnt always work 
#cat ${genome_prefix}*gff3 | grep "ID=gene" | grep -v "mRNA" > ${genome_prefix}.gene.gff
#echo cat ${genome_prefix}*gff3 | grep "ID=gene" | grep -v "mRNA" > ${genome_prefix}.gene.gff


# use bedtools to get the number of reads that map to specific regions

echo "bedtools count"
# for nomralisation use ? bedtools genomecov [OPTIONS] [-i|-ibam] -g (iff. -i)
cmd_count="bedtools multicov -bams ${genome_prefix}.bam -bed ${genome_prefix}.gene.gff > ${genome_prefix}_genomic.genes.cov"
echo ${cmd_count}
eval ${cmd_count}

cmd_ITS_count="bedtools multicov -bams ${genome_prefix}.bam -bed ${genome_prefix}.ITS.consensus.GFF > ${genome_prefix}_genomic.ITS.cov"
echo ${cmd_ITS_count}
eval ${cmd_ITS_count}

echo counting done


# just get the values using cut
# remove RNA genes from the all gene values, as these are what we are measuring in the other
# dataset. Assuming they are annotated in the GFF. It not... results may be a bit skewed.
cat ${genome_prefix}_genomic.genes.cov | grep -v "RNA" | cut -f10 > ${genome_prefix}_genomic.genes.cov.values
echo cat ${genome_prefix}_genomic.genes.cov | grep -v "RNA" | cut -f10 > ${genome_prefix}_genomic.genes.cov.values

cat  ${genome_prefix}_genomic.ITS.cov | cut -f10 >  ${genome_prefix}_genomic.ITS.cov.values
echo cat  ${genome_prefix}_genomic.ITS.cov | cut -f10 >  ${genome_prefix}_genomic.ITS.cov.values

# get stats summary of coverages
python ${path_to_ITS_clipping_file}/summary_stats.py --ITS ${genome_prefix}_genomic.ITS.cov.values --GFF ${genome_prefix}.ITS.consensus.GFF --all_genes_cov ${genome_prefix}_genomic.genes.cov.values -o ${genome_prefix}_stats.out
echo python ${path_to_ITS_clipping_file}/summary_stats.py --ITS ${genome_prefix}_genomic.ITS.cov.values --GFF ${genome_prefix}.ITS.consensus.GFF --all_genes_cov ${genome_prefix}_genomic.genes.cov.values -o ${genome_prefix}_stats.out



