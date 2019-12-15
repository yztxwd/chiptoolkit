#!/usr/bin/python python3.6
#coding = utf-8

"""
This module aims to perform pearson correlation between several datasets loaded from HDF5 cache
"""

import pandas as pd
import numpy as np
from scipy import stats
import sys
sys.path.append('/data1/yztxwd/Scripts/python/Module')
import numpyArrayDict as nad

from optparse import OptionParser, IndentedHelpFormatter

class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """loadHDF5_pearson.py

Load designated sources across the whole genome or in provided regions, then compute pearson/spearman coefficient. 

Example:
    $ python loadHDF5_pearson.py --specie mouse --database cache.h5 --source zero,one 
                                --region tss.region --bin 1000 -m pearson --prefix test [options]

Format Requirements:
    the first 5 columns of region file (1-based):
        chrom   start   end strand  id

@Copyright 2019, Jianyu Yang, Southern Medical University
"""

parser = OptionParser(description=description, formatter = CustomHelpFormatter())
parser.add_option('-c', '--specie', dest='specie', help='specie of data source')
parser.add_option('-d', '--database', dest='hdf5', help='hdf5 cache file')
parser.add_option('-s', '--source', dest='source', help='source index in hdf5 splitted by comma')
parser.add_option('-r', '--region', dest='region', default="", help='provide a region file to extract part of genome, default behavior is on whole genome')
parser.add_option('-b', '--bin', dest='bin', help='bin size to chop genomes into')
parser.add_option('-m', '--method', dest='method', help='pearson or spearman')
parser.add_option('-D', action='store_true', dest='dropZero', default=False, help='drop 0 count bins across all sources')
parser.add_option('-S', action='store_true', dest='storeBin', default=False, help='output file containing bin counts')
parser.add_option('-R', action='store_true', dest='reject', default=False, help='reject outliers')
parser.add_option('-p', '--prefix', dest='prefix', default='loadHDF5_pearson', help='output file prefix (suffix will be _$method.matrix)')
option, argument = parser.parse_args()

bin_size = int(option.bin)
sources = option.source.split(",")
num = len(sources)

# Return a boolean array with outliers marked with False
def reject_outliers(data, m = 1.5):
    Q3 = np.percentile(data, 75)
    Q1 = np.percentile(data, 25)
    IQR = Q3 - Q1
    idx = ~((data < (Q1 - 1.5 * IQR)) | (data > (Q3 + 1.5 * IQR)))
    return idx

# load region from hdf5 file, then concatenate them to array
def load_region(hdf5, region, source):
    returnArray = []
    region['keep'] = True
    for i in region.index:
        try:
            piece = hdf5.get_slice(source, region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end'])
        except:
            print("Exception when handling %s" %(region.loc[i, 'ID']))
            region.loc[i,'keep'] = False
            continue
        returnArray.extend(list(piece))
#    print("%s loaded" %(source))

    return returnArray

def do_pearson():
    # open hdf5 file
    hdf5 = nad.hdf5(option.hdf5, specie = option.specie)

    # Chop genomes into bins then concatenate and stack them
    bin_store = {}
    if option.region == "":
        for source in sources:
            nadict = hdf5.load_numpyArrayDict(source).store
            for i in nadict.keys():
                extra =  nadict[i].size%bin_size
                if(extra!=0):
                    append = bin_size - extra
                    nadict[i] = np.append(nadict[i], np.zeros(append))
                nadict[i] = nadict[i].reshape([-1, bin_size]).sum(axis=1)
            bin_store[source] = np.concatenate(list(nadict.values()))
    else:
        region = pd.read_csv(option.region, header=None, sep="\t", comment='#', usecols=[0,1,2,3,4])
        region = region.astype({1:'int',2:'int'})
        region.columns = ['chr', 'start', 'end', 'strand', 'ID']
        for source in sources:
            bin_store[source] = load_region(hdf5, region, source)
            extra = len(bin_store[source])%bin_size
            if(extra!=0):
                append = bin_size - extra
                bin_store[source].extend(np.zeros(append))
            bin_store[source] = np.array(bin_store[source]).reshape([-1, bin_size]).sum(axis=1)

    matrix = np.vstack(list(bin_store.values()))

    print(matrix.shape)
    
    # Drop 0 count bins in all sources if specified -D
    if(option.dropZero):
        idx = np.argwhere(np.all(matrix==0, axis=0))    
        matrix = np.delete(matrix, idx, axis=1)    

    # Reject outliers to avoid large value bias on correlation, run 10 times
    if(option.reject):
        for i in range(10):
            bArray = []
            for array in matrix:
                bArray.append(reject_outliers(array))
            bArray = np.vstack(bArray)
            matrix = matrix[:, np.any(bArray, axis=0)]

    # Compute pearson correlation matrix
    correlation_matrix = np.zeros(num**2).reshape(num, num)
    for i in range(num):
        for j in range(num):
            if option.method == 'pearson':
                correlation_matrix[i][j] = stats.pearsonr(matrix[i], matrix[j])[0]
            elif option.method == 'spearman':
                correlation_matrix[i][j] = stats.spearmanr(matrix[i], matrix[j])[0]

    # Retouch the matrix then output
    correlation_matrix = pd.DataFrame(correlation_matrix)
    correlation_matrix.index = correlation_matrix.columns = sources
    correlation_matrix.to_csv("%s_%s.matrix" %(option.prefix, option.method), header=True, sep="\t", index=True)
        
    # Store the bin counts if specified
    if(option.storeBin):
        bins = pd.DataFrame(bin_store)
        bins.to_csv("%s_bins.dataframe" %(option.prefix), header=True, sep='\t')
if __name__ == "__main__":
    do_pearson()


