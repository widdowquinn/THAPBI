#!/usr/bin/env python
#
# Tools for working with trimmomatic
#
# trimmomatic: http://www.bioinformatics.babraham.ac.uk/projects/trimmomatic/

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

    def run(self, L_reads, R_reads, outdir):
        """Run trimmomatic on the passed files"""
        self.__build_cmd(L_reads, R_reads, outdir)
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

    def __build_cmd(self, infname, outdir):
        """Build a command-line for trimmomatic"""
        self._outdirname = os.path.join(outdir, "trimmomatic_output")
        cmd = ["trimmomatic",
               infname,
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)
