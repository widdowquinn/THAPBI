#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to generate gff for ITS region BLAST hits)"
# The BLAST should already have been perfomed:
# blastn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out

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
        return file.read().split("\n")
    
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


def write_out_stats(ITS_cov, all_genes_cov, out):
    """function to write out summary stats. """
    # call function to get list of coverage per file.
    try:
        ITS_cov = parse_result_file(ITS_cov)
    except:
        raise ValueError("something wrong with ITS cov. file")
    try:
        all_genes_cov = parse_result_file(all_genes_cov)
    except:
        raise ValueError("something wrong with all genes cov. file")
    # out file to write to
    summary_stats_out = open(out, "w")
    title = "#min_cov\tmax_cov\tmean_cov\tstandard_dev\n"
    summary_stats_out.write(title)
    # call stats functio n
    ITSmin_cov, ITSmax_cov, ITSmean_cov, ITSstandard_dev = stat_tests(ITS_cov)
    ITS_data_formatted = "%s\tITSmin_cov%s\tITSmax_cov%s\tITSmean_cov%s\tITSstandard_dev\n" %(ITSmin_cov,\
                    ITSmax_cov, ITSmean_cov, ITSstandard_dev)
    summary_stats_out.write(ITS_data_formatted)

    GENEmin_cov, GENEmax_cov, GENEmean_cov, GENEstandard_dev = stat_tests(ITS_cov)
    GENE_data_formatted = "%s\tGENEmin_cov%s\tGENEmax_cov%s\tGENEmean_cov%s\tGENEstandard_dev\n" %(GENEmin_cov,\
                    GENEmax_cov, GENEmean_cov, GENEstandard_dev)
    summary_stats_out.write(GENE_data_formatted)

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

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", "--ITS", dest="ITS_cov", default=None,
                  help="the tab out file from the BLAST search",
                  metavar="FILE")
parser.add_option("-a", "--all_genes_cov", dest="all_genes_cov",
                  default=None",
                  help="the tab out file from the BLAST search",
                  metavar="FILE")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="stats.out",
                  help="outfile for the ITS and allgene stats")


(options, args) = parser.parse_args()


ITS_cov = options.ITS_cov
all_genes_cov = options.all_genes_cov
out_file = options.out_file



#run the program
file_list = ["ITS_cov"," all_genes_cov"]
for i in file_list:
    if not os.path.isfile(i):
        print("sorry, couldn't open the file: " + ex.strerror + "\n")
        print ("current working directory is :", os.getcwd() + "\n")
        print ("files are :", [f for f in os.listdir('.')])
        sys_exit("\n\nInput ITS file not found: %s" % i)

# call the top function    
write_out_stats(ITS_cov, all_genes_cov, out)


