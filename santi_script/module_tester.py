import errno
import logging
import logging.handlers
import multiprocessing
import os
import sys
import time
import traceback

from argparse import ArgumentParser

from thapbi_santi import trimmomatic

threads = 4

# Report last exception as string
def last_exception():
    """Returns last exception as a string, or use in logging."""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(traceback.format_exception(exc_type, exc_value,
                                              exc_traceback))


# Run as script
if __name__ == '__main__':


    # Set up logging
    logger = logging.getLogger('thapbi_santi_otus.py: %s' % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)

    try:
        logstream = open("module_testing.txt", 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     args.logfile)
        sys.exit(1)

    # Report input arguments
    logger.info("Command-line: %s" % ' '.join(sys.argv))
    logger.info("Starting testing: %s" % time.asctime())
    logger.info("starting trimmomatic testing")
    trim = trimmomatic.Trimmomatic("trimmomatic", logger)
    # this should break due to name
    #trimmomatic.trimmomatic("trimmomatic2", logger)


    logger.info("Trim reads by quality")
    trim.run("./data/DNAMIX_S95_L001_R1_001.fastq",
             "./data/DNAMIX_S95_L001_R2_001.fastq",
             threads, "trimmer_reads")
