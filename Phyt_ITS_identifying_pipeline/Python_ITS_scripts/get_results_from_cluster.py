#!/usr/bin/env python
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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import gzip



#####################################################################

def turn_fq_to_dic(in_fastq, barcode_length):
    """function to turn fq file into a dict
    read = seq[:length of bar code]
    return dict """
    read_to_barcode_dict = dict()
    #open the fastq file
    in_file = open(in_fastq)
    # iterate through the fatsq file
    try:
        for i, (title, seq, qual) in enumerate(FastqGeneralIterator(in_file)):
            read_to_barcode_dict[title.split(" ")[0]] = seq[:barcode_length]
    except:
        ValueError
        line_count = 0
        read_name = ""
        sequcene = ""
        for line in in_file:
            line_count = line_count+1
            if line.startswith("@"):
                #read line
                read_name = line
            if line_count % 4:
                #seq line
                sequcene = line[:barcode_length]
            read_to_barcode_dict[read_name.split(" ")[0]] = sequcene
    in_file.close()
    return read_to_barcode_dict


def get_barcode_seq (read, in_fastq):
    """function to get the barcode seq from the fq file.
    fq file needs to be unzipped.
    currently assuming barcode if 8bp and at the start of each seq??
    """
    return True


def parse_tab_file_get_clusters(filename1, left, right, barcode_length, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    #print coded_name_to_species_dict
    left_read_to_barcode_dic = turn_fq_to_dic(left, barcode_length)
    right_read_to_barcode_dic = turn_fq_to_dic(right, barcode_length)

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
        # are the database Phy singletons? - then we wont be interested
        number_of_reads_hitting_species = len(cluster_line)

        #set up some blank variables
        species = ""
        barcode_left = ""
        barcode_right = ""
        #database phy counter
        phy_count = 0
        for member in cluster_line:
            #is a memeber of the database in this cluster?
            if member.startswith("Phytophthora"):
                # yes we are interested in this cluster
                phy_count = phy_count+1
                #print member
                species = species+member+" "
                #remove all the memebr which are database memebers
                number_of_reads_hitting_species = number_of_reads_hitting_species-1
        if number_of_reads_hitting_species >1 and phy_count >0: #after the line, are there more memeber that are not database members?
            for member in cluster_line:
                if "Phytophthora" in member:
                    #no barcode for these - obviously!
                    continue
                print ("YYYEEEESSSSSS")
                #call func to get the bar codes
                left_barcode = left_read_to_barcode_dic[member]
                # add to str
                barcode_left = barcode_left+left_barcode+" "
                right_barcode = right_read_to_barcode_dic[member]
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

parser.add_option("-b", "--barcode_length", dest="barcode_length", default=8,
                  help="length of the barcode used. Default 8 ",
                  metavar="FILE")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

in_file = options.in_file
left = options.left
right = options.right
barcode_length = options.barcode_length
out_file = options.out_file

#run the program
barcode_length = int(barcode_length)

parse_tab_file_get_clusters(in_file, left, right, barcode_length, out_file)

print "done"
