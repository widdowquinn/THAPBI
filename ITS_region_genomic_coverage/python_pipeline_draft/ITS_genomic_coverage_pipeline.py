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
    #reads
    parser.add_argument("--left",, dest="left",
                        action="store", default=None,
                        help="left read file")
    parser.add_argument("--right",, dest="right",
                        action="store", default=None,
                        help="right read file")
    
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default=None,
                        help="Path to directory containing input reads")
    
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default="script_output",
                        help="Path to directory to write output")
    
    parser.add_argument("-r", "--reference", dest="reference_fasta",
                        action="store", default=None,
                        help="Path to reference sequence FASTA file")
    
    #programmes - can be in PATH
    parser.add_argument("--fastqc", dest="fastqc",
                        action="store", default="fastqc",
                        help="Path to FastQC executable")
    
    parser.add_argument("--blastn", dest="blastn",
                        action="store", default="blastn",
                        help="Path to blastn executable")
    
    parser.add_argument("--bowtie2", dest="bowtie2",
                        action="store", default="bowtie2",
                        help="Path to bowtie2 executable")
    
    parser.add_argument("--samtools", dest="samtools",
                        action="store", default="samtools",
                        help="Path to samtools executable")
    
    parser.add_argument("--bedtools", dest="bedtools",
                        action="store", default="bedtools",
                        help="Path to bedtools executable")
    
    #misc arguments
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

##########################################################################
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

    # Was a logfile specified? If so, use it
    if args.logfile is not None:
        try:
            logstream = open(args.logfile, 'w')
            err_handler_file = logging.StreamHandler(logstream)
            err_handler_file.setFormatter(err_formatter)
            err_handler_file.setLevel(logging.INFO)
            logger.addHandler(err_handler_file)
        except:
            logger.error("Could not open %s for logging" %
                         args.logfile)
            sys.exit(1)

    # Do we need verbosity?
    if args.verbose:
        err_handler.setLevel(logging.INFO)
    else:
        err_handler.setLevel(logging.WARNING)
    logger.addHandler(err_handler)

    # Report arguments, if verbose
    logger.info("Command-line: %s" % ' '.join(sys.argv))
    logger.info(args)
    logger.info("Starting pipeline: %s" % time.asctime())

    # Have we got an input directory, reference set and prefix? If not, exit.
        if args.indirname is None:
        logger.error("No input directory name (exiting)")
        sys.exit(1)
    logger.info("Input directory: %s" % args.indirname)
    if args.prefix is None:
        logger.error("No read file prefix given (exiting)")
        sys.exit(1)
    logger.info("Read file prefix: %s" % args.prefix)
    if args.reference_fasta is None:
        logger.error("No reference FASTA file given (exiting)")
        sys.exit(1)
    logger.info("Reference FASTA file: %s" % args.reference_fasta)

    # Have we got an output directory and prefix? If not, create it.
    if not os.path.exists(args.outdirname):
        logger.warning("Output directory %s does not exist - creating it" %
                       args.outdirname)
        os.makedirs(args.outdirname)
    # Check for the presence of space characters in any of the input filenames
    # If we have any, abort here and now.
    infilenames = sorted([os.path.join(args.indirname, fname) for
                          fname in os.listdir(args.indirname) if
                          fname.startswith(args.prefix)])
    logger.info("Input files: %s" % infilenames)
    for fname in infilenames:
        if ' ' in os.path.abspath(fname):
            logger.error("File or directory '%s' contains whitespace" % fname)
            logger.error("(exiting)")
            sys.exit(1)


##########################################################################
    # Check for presence of third-party tools, by instantiating interfaces
    logger.info("Checking third-party packages:")
    
    logger.info("\tFastQC... (%s)" % args.fastqc)
    fastQC = fastqc.FastQC(args.fastqc, logger)
    
    logger.info("\ttrimmomatic... (%s)" % args.trimmomatic)
    trimmomatic = trimmomatic.Trimmomatic(args.trimmomatic, logger)
    
    logger.info("\tblastn... (%s)" % args.blastn)
    blastn = blastn.Blastn(args.blastn, logger)
        
    logger.info("\tblastdb... (%s)" % args.blastn)
    makeblastdb = makeblastdb.MakeBlastdb(args.blastn, logger)

    logger.info("\tbowtie2... (%s)" % args.bowtie2)
    bowtie2 = bowtie2.Bowtie_map(args.bowtie2, logger)

    logger.info("\tsamtool... (%s)" % args.samtools)
    samtools_sort = samtools.Samtools_sort(args.samtools, logger)

    logger.info("\tsamtool... (%s)" % args.samtools)
    samtools_index = samtools.Samtools_index(args.samtools, logger)
    
    logger.info("\tbedtools... (%s)" % args.bedtools)
    betools = bedtools.Bedtools(args.bedtools, logger)

    # How many threads are we using?
    logger.info("Using %d threads/CPUs where available" % args.threads)
    
##########################################################################
    #now to run the commands.
    # Trim reads on quality - forward and reverse reads
    logger.info("Trim reads by quality")
    try:
        # this wont work. How to set it to produce the outfiles we require?
        trimmed_fnames = trimmomatic.run(args.left, args.right,\
                                         args.outdirname)
        logger.info("Trimmed FASTQ files:")
        logger.info("\t%s" % trimmed_fnames)
    except:
        logger.error("Error running trim_quality (exiting)")
        logger.error(last_exception())
        sys.exit(1)




    
