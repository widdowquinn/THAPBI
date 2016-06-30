#!/usr/bin/env python
#author: Peter Thorpe September 2016. The James Hutton Insitute,Dundee,UK.

#Title:
#script to reduce redundancy in GFF file)"

#imports
import os
import sys
from sys import stdin,argv
import sys
import datetime
from optparse import OptionParser


###########################################################################

def parse_tab_outfile(blast):
    """read in the blast tab file. Reads whole file into memory.
    returns a list, one list item per blast hit.
    """
    with open(blast) as file:
        return file.read().split("\n")

def split_gff_line(line):
    """function to split the gff line into its components"""
    scaf, genome, hit_number, start, stop, dot, direction, \
          dot2, description = line.split("\t")
    return scaf, genome, hit_number, start, stop, dot, direction, \
          dot2, description 
    


def write_out_ITS_GFF(gff, out): # this is a long function
    """function to write out the ITS blast hits in a GFF3
    like manner. """
    #this is out list of so called top/ longest matches which we will
    #append/remove as applicable
    blast_hit_to_info_dict = dict()
    best_stop = 0
    best_start = 0
    last_scaffold = "tmp"
    last_gff_line = "tmp"
    last_hit_number = ""
    
    # call function to get list of blast hits.
    try:
        blast_hits = parse_tab_outfile(gff)
    except:
        raise ValueError("something wrong with gff in file")
    GFF_out = open(out, "w")
    
    ############################################################################
    for i in blast_hits:
        #print "best_start, best_stop", best_start, best_stop
        if i.startswith("#"): #allows line to have comment.
            continue
        if not i.strip():
            continue #if the last line is blank
        
        assert len(i.split("\t"))== 9, ("""custom BLAST output?
                        not enough coloumns in gff file.""")
        #split up the gff line
        scaf, genome, hit_number, start, stop, dot, direction, \
                  dot2, description = split_gff_line(i)
        #populate the dictionaly
        blast_hit_to_info_dict[hit_number] = i
        
        # this is the first iteration. Populate the variables. 
        if last_scaffold == "tmp":
            best_stop = int(stop)
            best_start = int(start)
            last_scaffold = scaf
            last_gff_line = i
            last_hit_number = hit_number
            continue
        
        #lower and upper threshold for ranges. 
        lower_range_start = best_start - 50
        uppper_range_start = (best_stop - best_start) - 30
        
        upper_range_stop = best_stop + 50
        lower_range_stop = best_stop - best_start
        
        if scaf == last_scaffold:
            print "scaf", scaf, "last_scaffold", last_scaffold
            same_hit = False
            # same scaffold. Is the hit in the same region.
            #this means it could the same blast hit region
            if int(start) in range(lower_range_start, uppper_range_start):
                same_hit = True
                print "current :", i
                print "last time :", last_gff_line
                print "start :", start, "best_start ;", best_start, "\n\n"
                assert int(start) > int(best_start), ("""your
                    gff file has not been sorted by linux sort
                    please run this command
                    cat ${genome_prefix}.ITS.GFF | sort -k1,1 -k4,4 -k5,5
                    > sorted.gff""")
            #making sure it is the same blast hit region
            if int(stop) in range(lower_range_stop, upper_range_stop):
                same_hit = True

                # is the current stop greater or equal to the last?
                if int(stop) >= best_stop:
                    #print "stop" , stop, "best_stop", best_stop
                    #update the stop value
                    best_stop = stop
                    #remove this current dictionary entry
                    print ("iam removing: %s" %(hit_number))
                    same_hit = False # just to avoid removing again later
                    del blast_hit_to_info_dict[hit_number]

                    # get the old values again!
                    scaf, genome, hit_number, start, stop, dot, direction, \
                          dot2, description = split_gff_line(i)
                    #use this old vlaues with the new end coordinate
                    updated_values = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %(scaf,\
                                     genome, hit_number, \
                                     start, best_stop, \
                                     dot, direction, dot2, description)
                    
                    last_gff_line = updated_values
                    #update this blast enery in the 
                    blast_hit_to_info_dict[last_hit_number] = updated_values
            if same_hit:
                print ("iam removing: %s" %(hit_number))
                del blast_hit_to_info_dict[hit_number]
                
             # fill the variable at end of loop
            if not same_hit:
                best_stop = int(stop)
                best_start = int(start)
            last_scaffold = scaf
            last_gff_line = i
        else:
            best_stop = int(stop)
            best_start = int(start)
            last_scaffold = scaf
            last_gff_line = i
            last_hit_number = hit_number
    print_list = []        
    for key, val in blast_hit_to_info_dict.items():
        print_list.append(val)
    for i in sorted(print_list):
        print i
        
            
                    
                    
            
        
        
            
        
        

        #write to file
        #GFF_out.write(out_format)
    #close the write file
    #GFF_out.close()         


###########################################################################
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.1")
    sys.exit(0)


usage = """Use as follows:
Title:
script to reduce blast hits to consensus blast hit. This is required for
acurate genomic read coverage later on.

The BLAST should already have been perfomed:
 blastn -query ITS.fasta -db genome.fasta -outfmt 6 -out ITS_vs_geome.out
 a gff shoukld alread have been produced using 'generate_ITS_GFF.py'


$ filter_GFF.py --gff gff.out -o out.gff

Note:
coloumns 6,7 and 8 are 'made' up for the purpose of this task.
BLAST hits on the negative strand will be inverted so the
start is always less than the end coordinate.

"""

parser = OptionParser(usage=usage)

parser.add_option("-g", "--gff", dest="gff", default="outfmt6.out",
                  help="the tab out file from the BLAST search",
                  metavar="FILE")

parser.add_option("-o", "--out_file", dest="out_file",
                  default="ITS_GFF.out",
                  help="outfile for the ITS regions in GFF format")


(options, args) = parser.parse_args()


gff = options.gff
out_file = options.out_file



#run the program

if not os.path.isfile(gff):
    print("sorry, couldn't open the file: " + ex.strerror + "\n")
    print ("current working directory is :", os.getcwd() + "\n")
    print ("files are :", [f for f in os.listdir('.')])
    sys_exit("\n\nInput blast file not found: %s" % gff)

# call the top function    
write_out_ITS_GFF(gff, out_file)


