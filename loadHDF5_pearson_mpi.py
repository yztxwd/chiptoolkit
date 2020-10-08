#python3
#-*-coding:utf-8-*-

import sys
import numpyArrayDict as nad
import pandas as pd
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from mpi4py import MPI
from scipy import stats

from optparse import OptionParser, IndentedHelpFormatter

class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """loadHDF5_pearson_mpi.py

MPI version of loadHDF5_pearson.py
Load designated sources across the whole genome or in provided regions, then compute pearson/spearman coefficient. 

Example:
    $ mpiexec -n 4 python loadHDF5_pearson.py --specie mouse --database cache.h5 --source zero,one 
                                --region tss.region --bin 1000 -m pearson --prefix test [options]

Format Requirements:
    The first 5 columns of region file (1-based):
        chrom   start   end strand  id
    Blacklist should be bed file (0-based)

Dependencies:
    h5py-parallel (h5py built with --enable-parallel)
    mpi4py
    numpy
    pandas
    scipy

@Copyright 2019, Jianyu Yang, Southern Medical University
"""

parser = OptionParser(description=description, formatter = CustomHelpFormatter())
parser.add_option('-c', '--specie', dest='specie', help='specie of data source')
parser.add_option('-d', '--database', dest='hdf5', help='hdf5 cache file')
parser.add_option('-s', '--source', dest='source', help='source index in hdf5 splitted by comma')
parser.add_option('-r', '--region', dest='region', default="", help='provide a region file to extract part of genome, default behavior is on whole genome')
parser.add_option('-b', '--bin', dest='bin', help='bin size to chop genomes into')
parser.add_option('-m', '--method', dest='method', help='pearson or spearman')
parser.add_option('-B', '--blacklist', dest='blacklist', help='blacklist region')
parser.add_option('-D', action='store_true', dest='dropZero', default=False, help='drop 0 count bins across all sources')
parser.add_option('-S', action='store_true', dest='storeBin', default=False, help='output file containing bin counts')
parser.add_option('-R', action='store_true', dest='reject', default=False, help='reject outliers')
parser.add_option('-p', '--prefix', dest='prefix', default='loadHDF5_pearson', help='output file prefix (suffix will be _$method.matrix)')
option, argument = parser.parse_args()

bin_size = int(option.bin)
sources = option.source.split(",")
num = len(sources)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if num < size:
    raise Exception("Number of processors > Number of sources will cause unstability of program, please reduce core number")
    sys.exit(1)

# Global function definition
## Return a boolean array with outliers marked with True
def reject_outliers(data, m = 1.5):
    Q3 = np.percentile(data, 75)
    Q1 = np.percentile(data, 25)
    IQR = Q3 - Q1
    idx = ((data < (Q1 - m * IQR)) | (data > (Q3 + m * IQR)))
    return idx

## load region from hdf5 file, then concatenate them to array
def load_region(hdf5, region, source, blacklist):
    returnArray = []
    for i in region.index:
        try:
            if (not blacklist is None and sum(blacklist.get_range(region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end']))==0) or (blacklist is None):
                piece = hdf5.get_slice(source, region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end'])
                returnArray.extend(list(piece))
        except:
            print("Exception when handling %s" %(region.loc[i, 'ID']))
            continue
#    print("%s loaded" %(source))

    return returnArray

## load coverage sources in the given region, core function for workers
def do_pearson(hdf5, sources, region=None, blacklist=None):
    # Generate nad object for blacklist regions if supplied
    if not blacklist is None:
        blacklist = nad.numpyArrayDict(specie=option.specie).create_dict_fromDF(blacklist)

    # Chop genomes into bins then concatenate and stack them
    bin_store = {}
    if region is None:
        blacklistdict = blacklist.store
        for source in sources:
            nadict = hdf5.load_numpyArrayDict(source).store
            for i in nadict.keys():
                # Keep bases out of blacklist region
                nadict[i] = nadict[i][blacklistdict[i]==0]
                extra =  nadict[i].size%bin_size
                if(extra!=0):
                    append = bin_size - extra
                    nadict[i] = np.append(nadict[i], np.zeros(append))
                nadict[i] = nadict[i].reshape([-1, bin_size]).sum(axis=1)
            bin_store[source] = np.concatenate(list(nadict.values()))
    else:
        for source in sources:
            bin_store[source] = load_region(hdf5, region, source, blacklist)
            extra = len(bin_store[source])%bin_size
            if(extra!=0):
                append = bin_size - extra
                bin_store[source].extend(np.zeros(append))
            bin_store[source] = np.array(bin_store[source]).reshape([-1, bin_size]).sum(axis=1)
    matrix = np.vstack(list(bin_store.values()))
    matrix = matrix.astype('float64')
    return matrix

## Plot heatmap of correlation matrix
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """
    plt.figure(figsize=(10,10))

    # Plot heatmap
    df = pd.DataFrame(data)
    df.index = df.columns = sources
    g = sns.clustermap(df, cmap='viridis', linewidth=1)

    # Adjust plot
    bottom, top = g.ax_heatmap.get_ylim()
    g.ax_heatmap.set_ylim(bottom + 0.5, top - 0.5)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', rotation_mode="anchor")

# Load region in main processor, broadcast to workers
if rank == 0:
    if option.region:
        region = pd.read_csv(option.region, header=None, sep="\t", comment='#', usecols=[0,1,2,3,4])
        region = region.astype({1:'int',2:'int'})
        region.columns = ['chr', 'start', 'end', 'strand', 'ID']
    else:
        region = None
else:
    region = None
region = comm.bcast(region, root=0)

# Load blacklist regions if supplied, broadcast to workers
if rank == 0:
    blacklist = None
    if option.blacklist:
        blacklist = pd.read_csv(option.blacklist, header=None, sep='\t', comment='#', names=['chr', 'start', 'end'], usecols=[0,1,2])
        blacklist['start'] += 1
        blacklist['depth'] = 1
else:
    blacklist = None
blacklist = comm.bcast(blacklist, root=0)

# Assign jobs to each worker
if rank == 0:
    sources = option.source.split(',')
    interval = np.linspace(0, len(sources), size+1, dtype='int')
else:
    sources = None
    interval = None
interval = comm.bcast(interval, root=0)
sources = comm.bcast(sources, root=0)

# Time to work, send matrix to main processor
## Open HDF5 object for each processor
hdf5 = nad.hdf5(option.hdf5, specie = option.specie, mpi=True, mpi_comm=MPI.COMM_WORLD)
## Process each source
matrix = do_pearson(hdf5, sources[interval[rank]:interval[rank+1]], region = region, blacklist = blacklist)
collen = matrix.shape[1]
## Determine gather parameters
if rank == 0:
    recvbuf = np.empty([len(sources), collen], dtype='float64')
    split_sizes = []
    for i in range(len(interval)-1):
        split_sizes = np.append(split_sizes, interval[i+1]-interval[i])
    split_sizes_output = split_sizes * collen
    displacement_output = np.insert(np.cumsum(split_sizes_output), 0, 0)[0:-1]
else:
    split_sizes_output = None
    displacement_output = None
    recvbuf = None
split_sizes_output = comm.bcast(split_sizes_output, root=0)
displacement_output = comm.bcast(displacement_output, root=0)
## Send array to main processor
comm.Gatherv(sendbuf=matrix, recvbuf=[recvbuf, split_sizes_output, displacement_output, MPI.DOUBLE], root=0)
hdf5.close()
## Compute correlation, write to output
if rank == 0:   
    # Drop 0 count bins in all sources if specified -D
    if(option.dropZero):
        idx = np.argwhere(np.all(recvbuf==0, axis=0))    
        recvbuf = np.delete(recvbuf, idx, axis=1)    

    # Reject outliers to avoid large value bias on correlation
    if(option.reject):
        bArray = []
        for array in recvbuf:
            bArray.append(reject_outliers(array))
        bArray = np.vstack(bArray)
        recvbuf = recvbuf[:, ~np.any(bArray, axis=0)]

    # Compute correlation matrix
    correlation_matrix = np.zeros(num**2).reshape(num, num)
    for i in range(num):
        for j in range(num):
            try:
                if option.method == 'pearson':
                    correlation_matrix[i][j] = stats.pearsonr(recvbuf[i], recvbuf[j])[0]
                elif option.method == 'spearman':
                    correlation_matrix[i][j] = stats.spearmanr(recvbuf[i], recvbuf[j])[0]
            except ValueError:
                print(i, j)
                sys.exit(1)

    # Plot correlation heatmap
    heatmap(correlation_matrix, sources, sources)
    plt.savefig("%s_%s_heatmap.png" %(option.prefix, option.method), bbox_inches='tight', pad_inches=1, dpi=500)    

    # Retouch the matrix then output
    correlation_matrix = pd.DataFrame(correlation_matrix)
    correlation_matrix.index = correlation_matrix.columns = sources
    correlation_matrix.to_csv("%s_%s.matrix" %(option.prefix, option.method), header=True, sep="\t", index=True)
        
    # Store the bin counts if specified
    if(option.storeBin):
        bins = pd.DataFrame(recvbuf.T)
        bins.columns = sources
        bins.to_csv("%s_bins.dataframe" %(option.prefix), header=True, sep='\t')
