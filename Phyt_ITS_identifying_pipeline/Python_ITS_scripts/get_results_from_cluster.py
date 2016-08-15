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

def parse_line(line):
    """finction to parse line"""
    if not line.strip():
        return False #if the last line is blank
    if line.startswith("#"):
        return False
    #if not "Phytophthora" or "P." in line:
        #continue
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
    return cluster_line, number_of_reads_hitting_species, species, barcode_left, barcode_right, phy_count


def parse_tab_file_get_clusters(filename1, left, right, barcode_length, show_me_the_reads, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    #print coded_name_to_species_dict
    left_read_to_barcode_dic = turn_fq_to_dic(left, barcode_length)
    right_read_to_barcode_dic = turn_fq_to_dic(right, barcode_length)

    cluster_file = open (filename1, "r")
    summary_out_file = open(out_file, "w")
    #title for the results file
    if show_me_the_reads:
        title = "#species\tnumber_of_reads_hitting_species\tbarcodeR1\tbarcodeR2\treads_that_hit_species\n"
    else:
        title = "#species\tnumber_of_reads_hitting_species\tbarcodeR1\tbarcodeR2\n"
    summary_out_file.write(title)
    count = int(0)
    for line in cluster_file:
        count +=1
        #call function to parse line
        if not parse_line(line):
            continue
        reads_of_interest = ""
        cluster_line, number_of_reads_hitting_species, \
                      species, barcode_left, barcode_right, phy_count = parse_line(line)

        for member in cluster_line:
            member = member.rstrip()
            #is a memeber of the database in this cluster?
            if member.startswith("Phytophthora") or member.startswith("P."):
                # yes we are interested in this cluster
                phy_count = phy_count+1
                #print member
                species = species+member+" "
                #remove all the memebr which are database memebers
                number_of_reads_hitting_species = number_of_reads_hitting_species-1
        if number_of_reads_hitting_species >=1 and phy_count >0: #after the line, are there more memeber that are not database members?
            for member in cluster_line:
                if "Phytophthora" in member or "P." in member:
                    #no barcode for these - obviously!
                    continue
                print ("YYYEEEESSSSSS")
                print cluster_line
                # we do not add the reads that hit the 
                reads_of_interest = reads_of_interest+member+" "
                #call func to get the bar codes
                print "I am getting the left bar code", left_read_to_barcode_dic[member]
                try:
                    left_barcode = left_read_to_barcode_dic[member]
                except:
                    ValueError
                    left_barcode = "Not_found"
                # add to str
                barcode_left = barcode_left+left_barcode+" "
                print "I am getting the right bar code", right_read_to_barcode_dic[member]
                try:
                    right_barcode = right_read_to_barcode_dic[member]
                except:
                    ValueError
                    right_barcode = "Not_found"
                #add to str
                barcode_right = barcode_right+right_barcode+" "
            #format the data
            if show_me_the_reads:
                data_output = "%s\t%d\t%s\t%s\t%s\n" %(species.rstrip(), number_of_reads_hitting_species,\
                                                       barcode_left,\
                                                       barcode_right, reads_of_interest)
            else:
                data_output = "%s\t%d\t%s\t%s\t\n" %(species.rstrip(), number_of_reads_hitting_species,\
                                                       barcode_left,\
                                                       barcode_right)
            #write out the data
            summary_out_file.write(data_output)



    cluster_file.close()
    summary_out_file.close()
    return True

#############################################################################

#to run the script       

usage = """usage :

this scripts to summarise what clusters with what Phy species.


python parse_clusters.py -i clustering_outfile_already_decoded_from_temp_names -o summarise_clusters.out

command line option


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

parser.add_option("-s", "--show_me_the_reads", dest="show_me_the_reads", default=False,
                  help="show_me_the_reads in the output file for those that hit the Phy species."
                  " by default this is off, as the file could get very large. -s True if you want ...",
                  metavar="FILE")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()
#-i
in_file = options.in_file
# -l
left = options.left
# -r
right = options.right
# -b
barcode_length = options.barcode_length
# -o
out_file = options.out_file
# -s
show_me_the_reads = options.show_me_the_reads

#run the program
barcode_length = int(barcode_length)

parse_tab_file_get_clusters(in_file, left, right, barcode_length, show_me_the_reads, out_file)

print "done"
