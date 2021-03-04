import sys, os, csv
import numpy as np
import pandas as pd
from optparse import OptionParser, IndentedHelpFormatter
sys.path.append('/data2/yztxwd/scripts/python/Module')
import numpyArrayDict as nad


class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """loadHDF5.py

Load the coverage of sources in the given regions

Example:
    $ python loadHDF5.py --specie mouse --database cache.h5 --source H3K27me3 
                                    --region tss.region --prefix test -o ./

Format Requirements:
    The first 5 columns of region file (1-based):
        chrom   start   end strand  id

Dependencies:
    numpy
    pandas

@Copyright 2019, Jianyu Yang, Southern Medical University
"""

parser  = OptionParser(description= description, formatter = CustomHelpFormatter())
parser.add_option('-f','--file',dest='file',help='regions used to extract slices from hdf5')
parser.add_option('-c','--specie',dest='specie',help='species')
parser.add_option('-d','--database',dest='hdf5',help='hdf5 saving midpoint/coverage information')
parser.add_option('-s','--source',dest='source',help='source index in hdf5')
parser.add_option('-p', '--prefix',dest='prefix',default='loadHDF5',help='output file prefix')
parser.add_option('-o','--outDir',dest='outDir',default='.',help='output directory')
option, argument = parser.parse_args()

# get current working directory
cwd = os.getcwd()+"/"

# make directory
dirname = option.outDir
if (not os.path.exists(dirname)):
    os.makedirs(dirname)

# open hdf5 file
hdf5 = nad.hdf5(option.hdf5, specie=option.specie)

# load region file into pandas dataframe
region = pd.read_csv(option.file, header=None, sep="\t", comment='#', usecols=[0,1,2,3,4])
region = region[[0,1,2,3,4]]
region = region.astype({1:'int',2:'int'})
region.columns = ['chr', 'start', 'end', 'strand', 'ID']

# ready to extract region and write to output
ofile = open(dirname+'/'+option.prefix+'.csv', 'w')
writer = csv.writer(ofile, lineterminator='\n', delimiter='\t')
error = []
isFirst = True

writer.writerow(["#Source: %s" %(option.source)])
writer.writerow(["#Region: %s" %(option.file)])
for i in region.index:
    try:
        piece = hdf5.get_slice(option.source, region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end'])
        if region.loc[i,'strand'] == "-":
            piece = piece[::-1]
        ofile.write("%s\t" %(region.loc[i, 'ID']))
        writer.writerow(piece)
        if isFirst:
            sumRegion = piece
            isFirst = False
        else:
            sumRegion += piece
    except:
        error.append(region.loc[i, 'ID'])
ofile.close()
hdf5.close()
print("Error ID:\n%s" %(str(error)))

with open(dirname+'/'+option.prefix+'_sum.csv', 'w') as f:
    writer = csv.writer(f, lineterminator='\n', delimiter='\t')
    writer.writerow(["#Source: %s" %(option.source)])
    writer.writerow(["#Region: %s" %(option.file)])
    writer.writerow(sumRegion)





