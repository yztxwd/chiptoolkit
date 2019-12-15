#!/usr/bin/python3.6

import sys, os
from optparse import OptionParser, IndentedHelpFormatter
sys.path.append('/data2/yztxwd/scripts/python/Module')
import numpyArrayDict as nad


class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """
"""

parser  = OptionParser(description= description, formatter = CustomHelpFormatter())
parser.add_option('-f','--file',dest='file',help='midpoint/coverage bedgraph file separated by comma')
parser.add_option('-o','--specie',dest='specie',help='species')
parser.add_option('-d','--database',dest='hdf5',help='hdf5 saving midpoint/coverage information')
parser.add_option('-s','--source',dest='source',help='index name in hdf5 file')
option, argument = parser.parse_args()

filenames = option.file.split(',')

ifFirst = True
count = 0
for i in filenames:
    count += 1
    print("Processing No.%s file.........." %count)
    if ifFirst:
        intermediate = nad.numpyArrayDict(specie=option.specie).create_dict_fromfile(i)
        ifFirst=False
    else:
        intermediate.add_dict_fromfile(i)

print("Dumping into hdf5 file........")
intermediate.dump_hdf5(option.hdf5, option.source)
print("Script Done!")