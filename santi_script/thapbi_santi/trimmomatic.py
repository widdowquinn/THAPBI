#!/usr/bin/env python
#
# Tools for working with trimmomatic
#
# trimmomatic:
# http://www.bioinformatics.babraham.ac.uk/projects/trimmomatic/
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

# WARNING: This has not been tested at all YET. please do not use

import os
import sys

from subprocess import Popen, PIPE
from tools import is_exe


class Trimmomatic(object):
    """Class for working with trimmomatic"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            # TO DO: return this error message when it fails
            msg = ''.join(["trimmomatic executable not valid",
                   "trimming of reads will not be run...",
                   " SOLUTION: put trimmomatic in your PATH",
                   " and this is still failing "
                   "you may need to rename/ copy your trimmomatic ",
                   "binary to a file called timmomatic. ",
                   " please make sure this file is executable and",
                   "in your PATH",
                   " you can download the binaries from",
                   "http://www.usadellab.org/cms/?page=trimmomatic"])
            for m in msg:
                self._logger.warning(m)
            self._no_run = True
        self._exe_path = exe_path

# cmd_trimming="java -jar trimmomatic.jar PE -threads %d -phred33 ${left_read_file} ${right_read_file}
#${Name_of_project}_R1.fq.gz unpaired_R1.fq.gz ${Name_of_project}_R2.fq.gz unpaired_R2.fq.gz
# ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75" % (threads

    def __build_cmd(self, L_reads, R_reads, threads, outdir):
        """Build a command-line for trimmomatic"""
        self._outdirname = os.path.join(outdir, outdir)
        prefix = self._outdirname + "/" + L_reads.split("_R")[0]
        prefix = prefix.replace("./data", "")
        cmd = ["trimmomatic",
               "PE",
               "-threads", str(threads),
               "-phred33",
               L_reads, R_reads,
               prefix + "_paired_R1.fq.gz",
               "unpaired_R1.fq.gz",
               prefix + "_paired_R2.fq.gz",
               "unpaired_R2.fq.gz",
               "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10", "LEADING:3",
               "TRAILING:3", "SLIDINGWINDOW:4:25", "MINLEN:75"]
        self._cmd = ' '.join(cmd)

    def run(self, L_reads, R_reads, threads, outdir):
        """Run trimmomatic on the passed read files"""
        assert L_reads != R_reads, """Oh no,
        I am trying to perform trimming on two of the same files!
        Something has gone wrong in determining the left and right
        read files."""
        self.__build_cmd(L_reads, R_reads, threads, outdir)
        if not os.path.exists(self._outdirname):
            self._logger.info("Creating output directory: %s" %
                              self._outdirname)
            os.makedirs(self._outdirname)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        pipe = Popen(self._cmd, shell=True, stdout=PIPE)
        if pipe.wait() != 0:
            self._logger.error("trimmomatic generated some errors")
            sys.exit(1)
        return (self._outdirname, pipe.stdout.readlines())
