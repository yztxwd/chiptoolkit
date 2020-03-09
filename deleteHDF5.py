#!/usr/bin/python3.6

import sys, os, time
from optparse import OptionParser, IndentedHelpFormatter
sys.path.append("/data2/yztxwd/scripts/python/Module")
import numpyArrayDict as nad

class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """
This script aims to delete specific source in hdf5 file
"""

parser  = OptionParser(description= description, formatter = CustomHelpFormatter())
parser.add_option('-o','--specie',dest='specie',help='species')
parser.add_option('-d','--database',dest='hdf5',help='hdf5 saving midpoint/coverage information')
parser.add_option('-s','--source',dest='source',help='index name in hdf5 file')
option, argument = parser.parse_args()

# open hdf5 file
hdf5 = nad.hdf5(option.hdf5, specie = option.specie)

# delete specific source
print("""
Warning: %s will be deleted from %s
There will be 5 seconds before deletion, please be cautious about what you did
""" %(option.source, option.hdf5))
time.sleep(5)
hdf5.delete_source(option.source)

print("%s has been removed from %s" %(option.source, option.hdf5))
