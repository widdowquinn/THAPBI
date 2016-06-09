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
    parser.add_argument("--fastqc", dest="fastqc",
                        action="store", default="fastqc",
                        help="Path to FastQC executable")
    parser.add_argument("--trim_quality", dest="trim_quality",
                        action="store", default="trim_quality",
                        help="Path to seq_crumbs trim_quality script")
    parser.add_argument("--join_paired_ends", dest="join_paired_ends",
                        action="store", default="join_paired_ends.py",
                        help="Path to ea-utils join_paired_ends.py script")
    parser.add_argument("--convert_format", dest="convert_format",
                        action="store", default="convert_format",
                        help="Path to seq_crumbs convert_format script")
    parser.add_argument("--blastclust", dest="blastclust",
                        action="store", default="blastclust",
                        help="Path to blastclust")
    parser.add_argument("--muscle", dest="muscle",
                        action="store", default="muscle",
                        help="Path to MUSCLE")
    parser.add_argument("--pick_otus", dest="pick_otus",
                        action="store", default="pick_otus.py",
                        help="Path to QIIME pick_otus.py script")
    parser.add_argument("--pick_closed_reference_otus",
                        dest="pick_closed_reference_otus",
                        action="store",
                        default="pick_closed_reference_otus.py",
                        help="Path to QIIME pick_closed_reference_otus.py " +
                        "script")
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
