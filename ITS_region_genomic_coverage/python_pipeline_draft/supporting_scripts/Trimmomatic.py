#!/usr/bin/env python
#
# Tools for working with trimmomatic
#
# FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

import os
import sys

from subprocess import Popen, PIPE
from tools import is_exe


class trimmomatic(object):
    """Class for working with trimmomatic"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            msg = ["trimmomatic executable not valid",
                   "trimming of reads will not be run"]
            for m in msg:
                self._logger.warning(m)
            self._no_run = True
        self._exe_path = exe_path

    def run(self, infname_r1, infname_r2, outdir):
        """Run trimmomatic on the passed files"""
        self.__build_cmd(infname_r1, infname_r2, outdir)
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

    def __build_cmd(self, infname_r1, infname_r2, threads, outdir):
        """Build a command-line for trimmomatic"""
        self._outdirname = os.path.join(outdir, "trimmomatic_output")
        cmd =["java -jar trimmomatic-0.32.jar PE -threads",
        str(threads), "-phred33", infname_r1, infname_r2,
        "R1.fq.gz unpaired_R1.fq.gz R2.fq.gz unpaired_R2.fq.gz",
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3",
        "HEADCROP:9 TRAILING:3 SLIDINGWINDOW:4:22 MINLEN:51"]
        self._cmd = ' '.join(cmd)
