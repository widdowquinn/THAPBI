#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script perform stats on the coverage files already generated"

#imports
import os
import sys
import numpy
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################
# functions 
###########################################################################

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an interator
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

assert mean([1,2,3,4,5]) == 3

def parse_result_file(blast):
    """read in the blast tab file. Reads whole file into memeroy.
    returns a list, one list item per blast hit.
    """
    with open(blast) as file:
        data= file.read().split("\n")
        data1 = [line.rstrip("\n") for line in (data)
             if line.strip() != ""]
        return data1

def convert_to_int(results):
    """function to convert list of string to list
    of intergers"""
    results = [int(i) for i in results if i != ""]
    return results
    
def stat_tests(in_list):
    """function to return stats on a given list.
    returns min_cov, max_cov, mean_cov, standard_dev
    """
    min_cov = min(in_list)
    max_cov = max(in_list)
    mean_cov = mean(in_list)
    standard_dev = numpy.std(in_list)
    assert min_cov <= mean_cov <= max_cov
    return min_cov, max_cov, mean_cov, standard_dev


def write_out_stats(ITS_cov, GFF, all_genes_cov, out_file):
    """function to write out summary stats. """
    # call function to get list of coverage per file.
    number_of_ITS_blast_hits = len(parse_result_file(GFF))
    try:
        ITS_cov_str = parse_result_file(ITS_cov)
        #print ITS_cov_str
        ITS_cov = convert_to_int(ITS_cov_str)
    except:
        raise ValueError("something wrong with ITS cov. file")
    try:
        all_genes_cov_str = parse_result_file(all_genes_cov)
        all_genes_cov = convert_to_int(all_genes_cov_str)
    except:
        raise ValueError("something wrong with all genes cov. file")
    # out file to write to
    summary_stats_out = open(out_file, "w")
    title = "#gene_class\tmin_cov\tmax_cov\tmean_cov\tstandard_dev\n"
    summary_stats_out.write(title)
    
    # call stats function
    ITSmin_cov, ITSmax_cov, ITSmean_cov, ITSstandard_dev = stat_tests(ITS_cov)
    ITS_data_formatted = "ITS:\t%s\t%s\t%s\t%s\n" %(ITSmin_cov,\
                    ITSmax_cov, ITSmean_cov, ITSstandard_dev)
    #write out ITS results
    summary_stats_out.write(ITS_data_formatted)

    GENEmin_cov, GENEmax_cov, GENEmean_cov, GENEstandard_dev = stat_tests(all_genes_cov)
    GENE_data_formatted = "allGenes:\t%s\t%s\t%s\t%s\n" %(GENEmin_cov,\
                    GENEmax_cov, GENEmean_cov, GENEstandard_dev)
    summary_stats_out.write(GENE_data_formatted)

    blast_hits_info = "\nITS_blast_hit = %s \n" %(number_of_ITS_blast_hits)
    ratio_info = "ITS to gene ratio = %.1f \n" %(float(ITSmean_cov) / GENEmean_cov)
    summary_stats_out.write(blast_hits_info)
    summary_stats_out.write(ratio_info)
    

    #close the write file
    summary_stats_out.close()         


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
Title:
script to generate aummary stats for ITS coverage and all genes coverage


$ summary_stats.py --ITS ITS.cov --all all_gene.cov -o summary.out

ITS GFF file needed to count the number of ITS blast hits 
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", "--ITS", dest="ITS_cov", default=None,
                  help="coverage file for ITS regions",
                  metavar="FILE")
parser.add_option("-g", "--GFF", dest="GFF", default=None,
                  help="ITS GFF file",
                  metavar="FILE")
parser.add_option("-a", "--all_genes_cov", dest="all_genes_cov",
                  default=None,
                  help="the coverage file for all genes",
                  metavar="FILE")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="stats.out",
                  help="outfile for the ITS and allgene stats")


(options, args) = parser.parse_args()


ITS_cov = options.ITS_cov
GFF = options.GFF
all_genes_cov = options.all_genes_cov
out_file = options.out_file


#run the program
file_list = [ITS_cov, GFF, all_genes_cov]
for i in file_list:
    if not os.path.isfile(i):
        print("sorry, couldn't open the file: ", "\n")
        print ("current working directory is :", os.getcwd() + "\n")
        print ("files are :", [f for f in os.listdir('.')])
        sys.exit("\n\nInput ITS file not found: %s" % i)

# call the top function    
write_out_stats(ITS_cov, GFF, all_genes_cov, out_file)


