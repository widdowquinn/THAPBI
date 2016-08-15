#!/usr/bin/env python
#os importsimport os
from sys import stdin,argv
import sys
from optparse import OptionParser

#Title:
#script to trim barcodes or primer seq.

#############################################################################
#functions
def reformat_fasta_name(fasta, barcode, left_primer, right_primer, ITS, out):
    """function to retun the fa file with barcodes and primers removed"""
    f= open(out, 'w')
    f_in = open(fasta, "r")
    ITS = str(ITS)
    barcode = int(barcode)
    for seq_record in SeqIO.parse(fasta, "fasta"):
        #convert it to a string
        seq_record.seq= str(seq_record.seq)
        # make sure upper case
        seq_record.seq = seq_record.seq.upper()
        #remove the left and right barcode
        seq_record.seq = seq_record.seq[barcode:len(seq_record.seq)-barcode]
        if left_primer and right_primer:
            #if given left and right primer. Slice these off
            seq_record.seq = seq_record.seq[left_primer:len(seq_record.seq)-right_primer]
        #index the start of the ITS seq, slice there. 
        seq_record.seq = seq_record.seq[seq_record.seq.index(ITS):]
        assert seq_record.seq[:6] =="CCACAC"
        SeqIO.write(seq_record, f, "fasta")



#################################################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

$ python re_format_fasta.py -f seq.fasta --barcode 6 --left --right

default barcode length is 6. Change by passing an option

--left is the length of the left primer  - 0 by default
--right is the length of the right primer = 0 by default

GGAAGGTGAAGTCGTAACAAGG  = Primer ITS6 fwd

start of ITS = CCACAC

the script indexes for this instead.


requires Biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-f","--fasta", dest="fasta", default=None,
                  help="fasta file to have names altered")

parser.add_option("-b", "--barcode", dest="barcode", default=6,
                  help="length of barcodes",
                  metavar="FILE")

parser.add_option("-l", "--left", dest="left", default=False,
                  help="length of left primer",
                  metavar="FILE")

parser.add_option("-r", "--right", dest="right", default=False,
                  help="length right primer",
                  metavar="FILE")
parser.add_option("--ITS", dest="ITS", default="CCACAC",
                  help="length right primer",
                  metavar="FILE")

parser.add_option("-o", "--out", dest="out", default=None,
                  help="output filenames")


(options, args) = parser.parse_args()

fasta = options.fasta
barcode = options.barcode
left = options.left
right = options.right
ITS = options.ITS
out = options.out

#run the program

#biopython imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
reformat_fasta_name(fasta, barcode, left, right, ITS, out)

