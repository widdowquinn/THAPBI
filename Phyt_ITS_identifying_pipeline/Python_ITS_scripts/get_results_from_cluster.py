#title: Parse clusters and find the phy spcies
#in the cluster
#author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from optparse import OptionParser
import datetime
import os
from sys import stdin,argv
from Bio.SeqIO.QualityIO import FastqGeneralIterator



#####################################################################


def get_barcode_seq (read, in_fastq):
    """function to get the barcode seq from the fq file.
    fq file needs to be unzipped.
    currently assuming barcode if 8bp and at the start of each seq??
    """
    #open the fastq file
    in_file = open(in_fastq)
    # iterate through the fatsq file
    for i, (title, seq, qual) in enumerate(FastqGeneralIterator(in_file)):
        if read in title:
            in_file.close()
            return seq[:8]


def parse_tab_file_get_clusters(filename1, left, right, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    #print coded_name_to_species_dict
    cluster_file = open (filename1, "r")
    summary_out_file = open(out_file, "w")
    #title for the results file
    title = "#species\tnumber_of_reads_hitting_species\tbarcodeR1\tbarcodeR2\n"
    summary_out_file.write(title)
    count = int(0)
    for line in cluster_file:
        count +=1
        if not line.strip():
            continue #if the last line is blank
        if line.startswith("#"):
            continue
        if not "Phytophthora" or "P." in line:
            continue
        line = line.rstrip()
        out_put_str = ""
        cluster_line = line.rstrip("\n").split("\t")
        number_of_reads_hitting_species = len(cluster_line)
        species = ""
        barcode_left = ""
        barcode_right = ""
        #print cluster_line
        phy_count = 0
        for member in cluster_line:
            if member.startswith("Phytophthora"):
                phy_count = phy_count+1
                #print member
                species = species+member+" "
                #remove all the memebr which are database memebers
                number_of_reads_hitting_species = number_of_reads_hitting_species-1
        if phy_count >1:
            for member in cluster_line:
                if "Phytophthora" in member:
                    print ("YYYEEEESSSSSS")
                    #call func to get the bar codes
                    left_barcode = get_barcode_seq(member,left)
                    # add to str
                    barcode_left = barcode_left+left_barcode+" "
                    right_barcode = get_barcode_seq(member,right)
                    #add to str
                    barcode_right = barcode_right+right_barcode+" "
            data_output = "%s\t%d\t%s\t%s\t\n" %(species.rstrip(), number_of_reads_hitting_species,\
                                                       barcode_left,\
                                                       barcode_right)
            summary_out_file.write(data_output)



    cluster_file.close()
    summary_out_file.close()
    return True

#############################################################################

#to run the script       

usage = """usage :

this scripts return a 'master file' summarising what is in the clusters.


python parse_clusters.py -i clustering_outfile_already_decoded_from_temp_names -o summarise_clusters.out

command line option

requires:
biopython

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")

parser.add_option("-l", "--left", dest="left", default="R1.fq",
                  help="left reads unzipped fq file. Need this to get the barcode",
                  metavar="FILE")
parser.add_option("-r", "--right", dest="right", default="R2.fq",
                  help="right reads unzipped fq file. Need this to get the barcode",
                  metavar="FILE")
parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

in_file = options.in_file
left = options.left
right = options.right
out_file = options.out_file

#run the program

parse_tab_file_get_clusters(in_file, left, right, out_file)

print "done"
