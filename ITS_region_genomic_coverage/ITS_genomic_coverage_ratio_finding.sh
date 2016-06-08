#!/bin/bash
# put out files in current directory
#$ -cwd
#$ -l hostname="n13*"
# Deliver notifications to the following address
# Send notifications when the job begins (b), ends (e), or is aborted (a)
#$ -m a e -M peter.thorpe@hutton.ac.uk
#Abort on any error,
set -e
#(and thus quit the script right away)

echo Running on $HOSTNAME
echo Current PATH is $PATH
#source ~/.bash_profile

cd /home/pt40963/scratch/tree_health/
export TMP=/home/pt40963/scratch/${USER}_${JOB_ID}

##################################################################################################################################################################
# THESE VARIABLE NEED TO BE FILLED IN BY USER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

genome_prefix=Phytophthora_kernoviae.GCA_000333075.1.31.

genome_fasta=ftp://ftp.ensemblgenomes.org/pub/protists/release-31/fasta/phytophthora_kernoviae/dna/Phytophthora_kernoviae.GCA_000333075.1.31.dna.genome.fa.gz

genome_GFF=ftp://ftp.ensemblgenomes.org/pub/protists/release-31/gff3/phytophthora_kernoviae/Phytophthora_kernoviae.GCA_000333075.1.31.gff3.gz

read_1_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR278/008/SRR2785298/SRR2785298_1.fastq.gz

read_2_link=ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR278/008/SRR2785298/SRR2785298_2.fastq.gz

trimmomatic_path=/home/pt40963/Downloads/Trimmomatic-0.32

SRA_prefix=SRR2785298

path_to_ITS_clipping_file=~/misc_python/THAPBI/ITS_region_genomic_coverage

num_threads=12

# NOTHING TO BE FILLED IN BY USER FROM HERE!!!!
##################################################################################################################################################################


mkdir ${genome_prefix}
cd ${genome_prefix}


#reads genome and GFF, gunzip 

# i get my reads from: http://www.ebi.ac.uk/ena/data/warehouse/search
# use the advance.  Species ENTER, Lib.layout=paired,  platform=illumina, lib.source=genomic > search
wget ${read_1_link}
wget ${read_2_link}
wget ${genome_fasta}
wget ${genome_GFF}

# EXAMPLE: Phytophthora_kernoviae.GCA_000333075.1.31.dna.genome.fa.gz => Phytophthora_kernoviae.GCA_000333075.1.31.
gunzip ${genome_prefix}*


# blast to get representative ITS regions.
makeblastdb -in ${genome_prefix}*.fa -dbtype nucl
blastn -query ${path_to_ITS_clipping_file}/P.infestnas_ITS.fasta -db ${genome_prefix}*.fa -outfmt 6 -out n.Pi_ITS_vs_${genome_prefix}.out

# prepare ITS gff. python:

python ${path_to_ITS_clipping_file}/generate_ITS_GFF.py --blast n.Pi_ITS_vs_${genome_prefix}.out --prefix ${genome_prefix} -o ${genome_prefix}.ITS.GFF


#quality trim the reads
# head crop 9 as all illumina data seems to have non-ramdon bases at frist 9-12 nt. Weird!
java -jar ${trimmomatic_path}/trimmomatic-0.32.jar PE -threads ${num_threads} -phred33 ${SRA_prefix}_1.fastq.gz ${SRA_prefix}_2.fastq.gz R1.fq.gz unpaired_R1.fq.gz R2.fq.gz R2_unpaired.fq.gz ILLUMINACLIP:${path_to_ITS_clipping_file}/TruSeq3-PE.fa:2:30:10 LEADING:3 HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:51

rm  unpaired_R*

#map to the genome
mkdir $TMP

# index genome
bowtie2-build -f ${genome_prefix}*.fa bowtie_index_files

# randomly assign mutliple mapping reads ...  http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-p
#OLD command: bowtie2 --very-sensitive --non-deterministic --seed 1 --no-mixed --no-unal -p 8 -x Pi -1 R1.fq.gz -2 R2.fq.gz -S P.nicotiana.sam

#pipe SAM output on stdout into samtools to make it into BAM
#TODO: Put the temp unsorted BAM file on the local hard drive scratch space under /mnt/scratch
bowtie2 --very-sensitive --non-deterministic --seed 1 --no-mixed --no-unal -p ${num_threads} -x bowtie_index_files -1 R1.fq.gz -2 R2.fq.gz | samtools view -S -b -o $TMP/tmp_unsorted.bam -

# convert to sorted bam.
#samtools view -@ 8 -S -b -o tmp_unsorted.bam P.nicotiana.sam

samtools sort -@ ${num_threads} $TMP/tmp_unsorted.bam ${genome_prefix}
samtools index ${genome_prefix}.bam

# Clean up, doing this last in case fails
rm $TMP/tmp_unsorted.bam
# WARNING - Blindly deleting folders based on an encironemnt variable is DANGEROUS!
# There I am NOT doing: r m -r f $ T M P
# rmdir will only work if the directory is emtpy (and it should be)
rmdir $TMP

# use bedtools to get the number of reads that map to specific regions

bedtools multicov -bams ${genome_prefix}.bam -bed ${genome_prefix}.gene.gff > ${genome_prefix}_genomic.genes.cov

bedtools multicov -bams ${genome_prefix}.bam -bed ${genome_prefix}.ITS.GFF > ${genome_prefix}_genomic.ITS.cov


# just get the values using cut

cat ${genome_prefix}_genomic.genes.cov | cut -f10 > ${genome_prefix}_genomic.genes.cov.values

cat  ${genome_prefix}_genomic.ITS.cov | cut -f10 >  ${genome_prefix}_genomic.ITS.cov.values



