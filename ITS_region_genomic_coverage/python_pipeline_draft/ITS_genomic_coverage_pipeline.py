#!/usr/bin/env python
#
# ITS genomic coverage python pipeline
#
# Script to identify ITS regions and return genomic coverage.
#

#
# (c) The James Hutton Institute 2016
# Author: Peter Thorpe

import logging
import logging.handlers
import multiprocessing
import os
import sys
import time
import traceback

from argparse import ArgumentParser


from supporting_scripts import Fastqc, Trimmomatic, Blastn, \
     Summary_stats, Generate_ITS_GFF, MakeBlastdb.py\
     Bowtie_build, Bowtie_map, Samtools_sort, Samtools_index

##########################################################################
# make a folder for temp files
working_dir = os.getcwd()
dest_dir = os.path.join(working_dir, 'temp')
try:
    os.makedirs(dest_dir)
except OSError:
    print ("folder already exists, I will write over what is in there!!")
##########################################################################

# Process command-line arguments
def parse_cmdline(args):
    """Parse command-line arguments"""
    parser = ArgumentParser(prog="ITS_genomic_coverage_pipeline.py")
    parser.add_argument("-p", "--prefix", dest="prefix",
                        action="store", default=None,
                        help="Paired readfiles prefix")
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default=None,
                        help="Path to directory containing input reads")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default="script_output",
                        help="Path to directory to write output")
    parser.add_argument("-r", "--reference", dest="reference_fasta",
                        action="store", default=None,
                        help="Path to reference sequence FASTA file")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Report verbose output")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("-t", "--threads", dest="threads",
                        action="store", type=int,
                        default=multiprocessing.cpu_count(),
                        help="Number of threads to use (default: all)")
    return parser.parse_args()


# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))



# Run as script
if __name__ == '__main__':

    # Parse command-line
    args = parse_cmdline(sys.argv)

    # Set up logging
    logger = logging.getLogger('ITS_genomic_coverage_pipeline.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
