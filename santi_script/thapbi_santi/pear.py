#!/usr/bin/env python
#
# PEAR * (assemble overlapping reads)
# https://github.com/xflouris/PEAR
# follow this link to get the binaries. 
# http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz 
#
# (c) The James Hutton Institute 2016
# Author: Leighton Pritchard and Peter Thorpe

# WARNING: this module is not yet tested AT ALL.

import os
import sys

from subprocess import call
from tools import is_exe


class Pear(object):
    """Class for working with PEAR"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("""No PEAR program in PATH (exiting)
        The default name in the PEAR download is NOT PEAR. It is, for example
        pear-0.9.5-bin-64 . We recomment you change directory into the PEAR/bin
        forlder and 'cp pear-0.9.5-bin-64 pear' . This will create a copy of
        pear in that name. Make sure the PATH to this bin directory is in
        your PATH.

        If you are having troubles installing PEAR, pre-built binaries can be
        found at
        http://sco.h-its.org/exelixis/web/software/pear/files/pear-0.9.10-bin-64.tar.gz
        """)
            sys.exit(1)
        self._exe_path = exe_path

#cmd_pear="pear -f ${Name_of_project}_R1.fq.gz -r ${Name_of_project}_R2.fq.gz -o ${working_directory_path}/${Name_of_project}_outfiles/${Name_of_project}_PEAR" 
    def run(self, L_reads, R_reads, outdir):
        """Run PEAR on the passed files"""
        self.__build_cmd(L_reads, R_reads, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("PEAR terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("pick_otus.py returned %s" % retcode)
        return self._outdirname

    def __build_cmd(self, L_reads, R_reads, outdir):
        """Build a command-line for pick_otus.py"""
        self._outdirname = os.path.join(outdir, "qiime_uclust_OTUs")
        cmd = ["pick_otus.py",
               "-m", "uclust_ref",
               "-s", "0.99",
               "-i", infname,
               "-r", refname,
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)


class Pick_Closed_Ref_Otus(object):
    """Class for working with pick_closed_reference_otus.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No pick_closed_reference_otus.py script " +
                               "(exiting)")
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, infnames, refname, outdir):
        """Run pick_closed_reference_otus.py on the passed file"""
        self.__build_cmd(infnames, refname, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("pick_closed_reference_otus.py terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("pick_closed_reference_otus.py returned " +
                              "%s" % retcode)
        return self._outdirname

    def __build_cmd(self, infname, refname, outdir):
        """Build a command-line for pick_closed_reference_otus.py"""
        self._outdirname = os.path.join(outdir, "qiime_closedref_OTUs")
        cmd = ["pick_closed_reference_otus.py", "-f",
               "-i", infname,
               "-r", refname,
               "-o", self._outdirname]
        self._cmd = ' '.join(cmd)


class Join_Paired_Ends(object):
    """Class for working with join_paired_ends.py"""
    def __init__(self, exe_path, logger):
        """Instantiate with location of executable"""
        self._logger = logger
        self._no_run = False
        if not is_exe(exe_path):
            self._logger.error("No join_paired_ends.py script (exiting)")
            sys.exit(1)
        self._exe_path = exe_path

    def run(self, infnames, outdir):
        """Run joined_paired_ends.py on the passed file"""
        self.__build_cmd(infnames, outdir)
        msg = ["Running...", "\t%s" % self._cmd]
        for m in msg:
            self._logger.info(m)
        retcode = call(self._cmd, shell=True)
        if retcode < 0:
            self._logger.error("join_paired_ends.py terminated by " +
                               "signal %s" % -retcode)
            sys.exit(1)
        else:
            self._logger.info("join_paired_ends.py returned %s" % retcode)
        return self._outfname

    def __build_cmd(self, infnames, outdir):
        """Build a command-line for join_paired_ends.py"""
        f1, f2 = tuple(infnames)
        cmd = ["join_paired_ends.py",
               "-f", f1,
               "-r", f2,
               "-o", outdir]
        self._cmd = ' '.join(cmd)
        self._outfname = os.path.join(outdir, "fastqjoin.join.fastq")
