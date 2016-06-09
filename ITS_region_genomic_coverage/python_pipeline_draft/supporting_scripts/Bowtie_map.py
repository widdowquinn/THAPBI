#!/usr/bin/env python
#
# Tools for working with BLASTn
#
# (c) The James Hutton Institute 2016
# Author: Peter Thorpe

import os
import sys

from subprocess import Popen, PIPE
from tools import is_exe


class Bowtie2_build(object):
    """Class for working with bowtie2-build"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No bowtie2-build at %s (exiting)" % exe_path)
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, infnames, outdir, threads):
        """Run bowtie2-build on the passed file"""
        self.__build_cmd(infnames, outdir, threads)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        pipe = Popen(self._cmd, shell=True, stdout=PIPE)
        if pipe.wait() != 0:
            self._logger.error("bowtie2-build generated some errors")
            sys.exit(1)
        return (self._outfname, pipe.stdout.readlines())

    def __build_cmd(self, infname, outdir, threads):
        """Build a command-line for bowtie2-build"""
        
        cmd = ["bowtie2 --very-sensitive --non-deterministic --seed 1",
               "--no-mixed --no-unal -p", str(threads),
               infname,
               "-x bowtie_index_files -1 R1.fq.gz -2 R2.fq.gz | samtools view -@",
               str(threads),"-S -b -o temp/tmp_unsorted.bam" ]
        self._cmd = ' '.join(cmd)
