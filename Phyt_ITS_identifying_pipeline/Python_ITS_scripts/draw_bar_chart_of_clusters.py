#title: Parse clusters and find the corresponding spcies
#in the cluster
#author: Peter Thorpe September 2016. The James Hutton Insitute, Dundee, UK.

#imports
import os
import sys
from optparse import OptionParser
import datetime
import os
from sys import stdin,argv
#imports for graphs
#import seaborn as sns
import matplotlib
import numpy


##############################################################################################
# drawing a histogram
#n, bins, patches = hist(data)

def parse_line(line):
    """function to parse a given line and return
    tab separated elements"""
    if not line.strip():
        return False #if the last line is blank
    if line.startswith("#"):
        return False
    cluster_line_split = line.rstrip("\n").split("\t")
    #print ("I assume the element are tab separated")
        #cluster_line_split = line.rstrip("\n").split()
    return cluster_line_split
    

def count_element_in_cluster(cluster_line_split):
    """func to count to total memebers in each cluster.
    Returns a dic with """
    species_set = set([])
    #counters to keep track
    members_count = 0
    species_count = 0
    for member in cluster_line_split:
        members_count = members_count+1
        if member not in species_set:
            species_count = species_count+1
            species_set.add(member)
    return members_count, species_count

try:
    # New in Python 3.4
    from statistics import mean
except ImportError:
    def mean(list_of_values):
        """Calculate the mean average of a list of numbers."""
        # Quick and dirty, assumes already a list not an interator
        # so don't have to worry about getting the divisor.
        # Explicit float(...) to allow for Python 2 division.
        return sum(list_of_values) / float(len(list_of_values))

assert mean([1,2,3,4,5]) == 3
def covert_dict_to_list_of_value(in_dict):
    """function to convert a given dict, to convert it to a
    list of value. In dict should be something like this:
    count_dict= {1: 4, 2: 1}  """
    #need this value to plot the number of histogram bars
    number_of_keys = len(in_dict.keys())
    output_list = []
    for key, val in in_dict.items():
	for i in range(0,val):
            output_list.append(key)
    return sorted(output_list), number_of_keys

def plot_graph(data_values, title, number_of_keys, file_in):
        #import matplotlib
        matplotlib.use('Agg')
        import pylab
        import matplotlib.mlab as mlab
        bins = max(data_values)
        pylab.hist(data_values, facecolor='blue')
        mu = mean(data_values)
        standard_dev = numpy.std(data_values)
        # add a 'best fit' line
        y = mlab.normpdf(bins, mu, standard_dev)
        l = pylab.plot(bins, y, 'r--', linewidth=10)
        pylab.grid(True)
        pylab.title(title)
        pylab.xlabel('number in cluster')
        pylab.ylabel('Count')
        pylab.show()
        pylab.savefig(file_in+"_"+title+'_histogram.png')

        os.chdir('..')

    
def parse_tab_file_get_clusters(filename1, out_file):
    """#script to open up a tab separeted clustering output and identify the
    species in the clustering"""
    cluster_file = open (filename1, "r")
    summary_out_file = open(out_file, "w")
    #dictionaries for keeping the counts
    member_in_cluster_to_count_dict = dict()
    species_in_cluster_count_dict = dict()

    # a way of keeping track of the iteration
    interation_count = int(0)
    #iterate through the file
    for line in cluster_file:
        interation_count +=1
        #call the func to split up the line
        cluster_line_split = parse_line(line.rstrip())
        if not cluster_line_split:
            #this could be a blank line or a line that starts with #
            continue
        
        # call the function to get the number of elements and spceis. 
        members_count, species_count = count_element_in_cluster(cluster_line_split)
        try:
            #if we have seen this count before, then just add one to it.
            member_in_cluster_to_count_dict[members_count] +=1
        except:
            KeyError
            #print ("not seen this before")
            #not seen this before, set up a new dic element and make the equal 1
            member_in_cluster_to_count_dict[members_count] = 1
        try:
            #if we have seen this count of species before,
            #then just add one to it.
            species_in_cluster_count_dict[species_count] +=1
        except:
            KeyError
            species_in_cluster_count_dict[species_count] = 1
        

    #print ("species_in_cluster_count_dict", species_in_cluster_count_dict)
    #print ("member_in_cluster_to_count_dict", member_in_cluster_to_count_dict)

    #call the function to convert dic to list
    species_in_cluster_list, species_number_of_keys = covert_dict_to_list_of_value(species_in_cluster_count_dict)
    member_in_cluster_list, member_number_of_keys = covert_dict_to_list_of_value(member_in_cluster_to_count_dict)
    #print ("we have this KEYS species in cluster", species_number_of_keys)
    #print ("we have this KEYS memebers in cluster", member_number_of_keys)

    plot_graph(species_in_cluster_list, "species_in_cluster", species_number_of_keys, filename1)
    plot_graph(member_in_cluster_list, "member_in_cluster", member_number_of_keys, filename1)


    return True

#############################################################################

#to run the script       

usage = """usage :

script to graphically represent clustering output.


python draws_bar...py -i clustering_file -o summarise_clusters.out

requires:
use pip install ...
 seaborn
 matplotlib
 numpy

"""

parser = OptionParser(usage=usage)

parser.add_option("-i","--in", dest="in_file", default=None,
                  help="clustering out file")

parser.add_option("-o", "--out_prefix", dest="out_file", default="summarise_clusters.out",
                  help="prefix to the output filenames")


(options, args) = parser.parse_args()

in_file = options.in_file
out_file = options.out_file

#run the program

parse_tab_file_get_clusters(in_file, out_file)

print "done"
