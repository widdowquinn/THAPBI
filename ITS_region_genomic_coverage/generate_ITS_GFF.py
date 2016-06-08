#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to generate gff for ITS region BLAST hits)"
# The BLAST should already have been perfomed:
# blastn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out

#imports
import os
import sys
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################

def parse_blast_tab_outfile(blast):
    """read in the blast tab file. Reads whole file into memeroy.
    returns a list, one list item per blast hit.
    """
    with open(blast) as file:
        return file.read().split("\n")
    
def spit_blast_data(i):
    """function to split up the blast hits
    and return \t formatted data. checks start < stop
    , alters these if need be..."""
    # split the blast line and assign the feilds respectively   
    queryId, subjectId, percIdentity, alnLength,mismatchCount,\
             gapOpenCount, queryStart, queryEnd, subjectStart, \
             subjectEnd, eVal, bitScore = i.split("\t")
    #reverse negative blast hits (breaks bamtools if not fixed)
    if int(subjectStart) > int(subjectEnd):
        temp_subjectStart = subjectEnd
        temp_subjectEnd = subjectStart
        out_format="%s\t%s\tITS_blast_hit_%d\t"
            "%s\t%s\t.\t+\t.\tITS_blast_hits_region" \
              %(subjectId, prefix, blast_count, \
                temp_subjectStart,temp_subjectEnd)
    else:
        #direction find. ready for writing out. 
        out_format= "%s\t%s\tITS_blast_hit_%d\t"
            "%s\t%s\t.\t+\t.\tITS_blast_hits_region"\
              %(subjectId, prefix, blast_count, subjectStart,\
                subjectEnd)
    return out_format


def write_out_ITS_GFF(blast, prefix, out):
    """function to write out the ITS blast hits in a GFF3
    like manner. """
    # call function to get list of blast hits.
    try:
        blast_hits = parse_blast_tab_outfile(blast)
    except:
        raise ValueError:
            print("something wrong with the blast out file")
    GFF_out = open(out, "w")
    # counter to index the blast hits in the GFF file
    blast_count = 0
    # santity check to remove duplicate events
    already_seen_set = set([])

    for i in blast_hits:
        if i.startswith("#"):
            #allows the outfile to have comment lines.
            continue
        blast_count = blast_count +1
        # check this is a unique blast hit. Remove duplicates!
        if i not in already_seen_set:
            #add this to seen set. 
            already_seen_set.add(i)
            if len(i.split("\t")) > 12:
                #remove tax id and extra coloumns - not needed.
                i = i[:12]
            if len(i.split("\t")) >12:
                raise ValueError:
                    print ("custom BLAST output? not enough "
                           "coloumns in blast file. ")
        out_format = spit_blast_data(i)
        #write to file
        GFF_out.write(out_format)
    #close the write file
    GFF_out.close()         


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print "v0.0.1"
    sys.exit(0)


usage = """Use as follows:
Title:
script to generate gff for ITS region BLAST hits)"
The BLAST should already have been perfomed:
 blastn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out


$ generate_ITS_GFF.py -b blast.out --prefix p.infestans -o gff.out

Note:
coloumns 6,7 and 8 are 'made' up for the purpose of this task.
BLAST hits on the negative strand will be inverted so the
start is always less than the end coordinate.

"""

parser = OptionParser(usage=usage)

parser.add_option("-b", "--blast", dest="blast", default="outfmt6.out",
                  help="the tab out file from the BLAST search",
                  metavar="FILE")
parser.add_option("--prefix", dest="prefix",
                  default="temp_name",
                  help="name for column 2 in GFF. Best to "
                  " use the origin of the data")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="ITS_GFF.out",
                  help="outfile for the ITS regions in GFF format")


(options, args) = parser.parse_args()


blast = options.blast
prefix = options.prefix
out_file = options.out_file



#run the program

if not os.path.isfile(blast):
    print("sorry, couldn't open the file: " + ex.strerror + "\n")
    print ("current working directory is :", os.getcwd() + "\n")
    print ("files are :", [f for f in os.listdir('.')])
    sys_exit("\n\nInput blast file not found: %s" % blast)

# call the top function    
write_out_ITS_GFF(blast, prefix, out_file)


