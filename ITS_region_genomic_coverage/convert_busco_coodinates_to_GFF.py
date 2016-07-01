#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to generate gff for ITS region busco_coord hits)"
# The busco_coord should already have been perfomed:
# busco_coordn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out

#imports
import os
import sys
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################

def parse_tab_outfile(busco_coord):
    """read in the busco_coord tab file. Reads whole file into memeroy.
    returns a list, one list item per busco_coord hit.
    """
    with open(busco_coord) as file:
        return file.read().split("\n")



def write_out_ITS_GFF(busco_coord, prefix, out):
    """function to write out the ITS busco_coord hits in a GFF3
    like manner. """
    # call function to get list of busco_coord hits.
    try:
        busco_coord_hits = parse_tab_outfile(busco_coord)
    except:
        raise ValueError("something wrong with busco_coord out file")
    GFF_out = open(out, "w")
    # counter to index the busco_coord hits in the GFF file
    busco_coord_count = 0
    # santity check to remove duplicate events
    already_seen_set = set([])

    for i in sorted(busco_coord_hits):
        if i.startswith("#"):
            #allows the outfile to have comment lines.
            continue
        if not i.strip():
            continue #if the last line is blank
        busco_coord_count = busco_coord_count +1
        # check this is a unique busco_coord hit. Remove duplicates!
        EOG, scaffold, start, stop = i.split("\t")
        
        out_format = "%s\t%s\tEOG_hit_%d\t%s\t%s\t.\t+\t.\t%s\n" %(scaffold,\
                            prefix, busco_coord_count,\
                            start, stop, EOG)
        #write to file
        GFF_out.write(out_format)
    #close the write file
    GFF_out.close()         


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
Title:
script to generate gff for  busco_coord hits"
The busco_coord should already have been perfomed:
 busco_coordn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out


$ convert_busco_coodinates_to_GFF.py -b busco_coord.out --prefix p.infestans -o gff.out

Note:
coloumns 6,7 and 8 are 'made' up for the purpose of this task.
busco_coord hits on the negative strand will be inverted so the
start is always less than the end coordinate.

"""

parser = OptionParser(usage=usage)

parser.add_option("-b", "--busco_coord", dest="busco_coord", default="outfmt6.out",
                  help="the tab out file from the busco_coord search",
                  metavar="FILE")
parser.add_option("--prefix", dest="prefix",
                  default="temp_name",
                  help="name for column 2 in GFF. Best to "
                  " use the origin of the data")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="ITS_GFF.out",
                  help="outfile for the ITS regions in GFF format")


(options, args) = parser.parse_args()


busco_coord = options.busco_coord
prefix = options.prefix
out_file = options.out_file



#run the program

if not os.path.isfile(busco_coord):
    print("sorry, couldn't open the file: " + ex.strerror + "\n")
    print ("current working directory is :", os.getcwd() + "\n")
    print ("files are :", [f for f in os.listdir('.')])
    sys_exit("\n\nInput busco_coord file not found: %s" % busco_coord)

# call the top function    
write_out_ITS_GFF(busco_coord, prefix, out_file)


