#!/usr/bin/env python
#os importsimport os
from sys import stdin,argv
import sys
from optparse import OptionParser

#Title:
#script to reformt fasta ID names.

#############################################################################
#functions
def reformat_fasta_name(filename, databse, out):
    "this function re-write a file as a fasta file but with altered names"
    f= open(out, 'w')
    f_in = open(fasta, "r")

    name = []
    count = 0
    for seq_record in SeqIO.parse(filename, "fasta"):
        count = count+1
        old_to_new_name = "%s\ttemp_seq_%d\n" % (seq_record.id, count)
        name.append(old_to_new_name)
        #remove the read prefix to get to uniq names
        seq_record.id = "temp_seq_%d" % (count)
        seq_record.description = ""
        SeqIO.write(seq_record, f, "fasta")
    name_out = open(databse, "w")
    #file to keep track of the original names if we need them
    for i in name:
       name_out.write(i)
    name_out.close()
    f.close()
    f_in.close()




#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python complete....py -f seq.fasta --prefix @M01157:20:000000000-D07KA:

script to reformt fasta names. The read name can bre

requires Biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="fasta file to have names altered")

parser.add_option("-d", "--databse", dest="databse", default=None,
                  help="outfile to keep track of old and new ID names",
                  metavar="FILE")
parser.add_option("-o", "--out", dest="out", default=None,
                  help="output filename")


(options, args) = parser.parse_args()

fasta = options.fasta
databse = options.databse
out = options.out

#run the program

#biopython imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
reformat_fasta_name(fasta, databse, out)

