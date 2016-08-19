
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

def name_set(in_names):
    name_set = set([])
    for i in in_names:
        if not i.startswith("#"):
            name_set.add(i)
    return name_set
    

def reformat_as_fasta(filename, outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')
    for seq_record in SeqIO.parse(filename, "fasta"):
        seq_record.seq = str(seq_record.seq)
        seq_record.seq = seq_record.seq.replace("-" , "")
        SeqIO.write(seq_record, f, "fasta")                    
    
    f.close()
    return True

def reformat_as_fasta_nobio(filename, outfile):
    "this function re-write a file as a fasta file"
    f= open(outfile, 'w')
    fh = open(filename, 'r')
    for line in fh:
        if line.startswith(">"):
            title = line.split("/")[0]
            print >> f, title.rstrip()
        else:
            seq = line.replace("-" , "")
            print >> f, seq.rstrip()
    
    f.close()
    return True



if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:

converts

$ python rewrite_as_fasta.py -i in.fasta -l min_length_of_seq (default(3)) --not_wanted --wanted -o out.fasta

script either reformats badly formated fasta file. Within reason. Ie. Word documents will still break it.

if lenght of seq is longer than -l default 3 - writes to file.

can filter fasta by giving it a list of --not_wanted or --wanted names.

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="current fasta you want to reformat")


parser.add_option("-o", "--out", dest="out", default=None,
                  help="Output filename",
                  metavar="FILE")




(options, args) = parser.parse_args()

in_file = options.in_file

out = options.out



#reformat_as_fasta(in_file,out)
reformat_as_fasta_nobio(in_file, out)
print 'done'

