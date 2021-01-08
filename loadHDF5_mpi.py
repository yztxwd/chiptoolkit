import sys, os, csv
import numpy as np
import pandas as pd
import h5py
import matplotlib
import matplotlib.pyplot as plt
from mpi4py import MPI
from optparse import OptionParser, IndentedHelpFormatter
from collections import defaultdict
#sys.path.append('/data2/yztxwd/scripts/python/Module')
import numpyArrayDict as nad


class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """loadHDF5_mpi.py

MPI version of loadHDF5.py
Load the coverage of sources in the given regions

Example:
    $ mpiexec -n 4 python loadHDF5_mpi.py --specie mouse --database cache.h5 --source zero,one
                                    --region tss.region --prefix test -o ./

Format Requirements:
    The first 5 columns of region file (1-based):
        chrom   start   end strand  id
    Blacklist should be bed file (0-based)

Dependencies:
    h5py-parallel (h5py build with --enable-parallel)
    mpi4py
    numpy
    pandas

@Copyright 2019, Jianyu Yang, Southern Medical University
"""

parser  = OptionParser(description= description, formatter = CustomHelpFormatter())
parser.add_option('-r','--region',dest='region',help='regions used to extract slices from hdf5')
parser.add_option('-c','--specie',dest='specie',help='species')
parser.add_option('-d','--database',dest='hdf5',help='hdf5 saving midpoint/coverage information')
parser.add_option('-s','--source',dest='source',help='source index in hdf5')
parser.add_option('-b','--blacklist',dest='blacklist',help='blacklist region')
parser.add_option('-n','--normalize',dest='normalize',default=None,help="Method for normalization")
parser.add_option('-p', '--prefix',dest='prefix',default='loadHDF5',help='output file prefix')
parser.add_option('-o','--outDir',dest='outDir',default='.',help='output directory')
option, argument = parser.parse_args()

sources = [s.strip() for s in option.source.split(",")]
num = len(sources)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# make directory
dirname = option.outDir
if rank == 0:
    if (not os.path.exists(dirname)):
        os.makedirs(dirname)

if num < size:
    raise Exception("Number of processors > Number of sources will cause unstability of program, please reduce core number")
    sys.exit(1)

# Global function definition
## set outlier to a smaller value
def shrink_outlier(data, m=3):
    Q3 = np.percentile(data, 75)
    Q1 = np.percentile(data, 25)
    IQR = Q3 - Q1
    idx1 = data < (Q1 - m * IQR)
    idx2 = data > (Q3 + m * IQR)
    data[idx1] = (Q1 - m * IQR)
    data[idx2] = (Q3 + m * IQR)
    return data

def toPercentile(a, dic):
    return np.vectorize(dic.__getitem__)(a)

def get_coverage(hdf5, sources, region, blacklist):
    # get current working directory
    cwd = os.getcwd()+"/"

    # generate nad object if blacklist supplied
    if not blacklist is None:
        print("Use blacklist to exclude region")
        blacklist = nad.numpyArrayDict(specie=option.specie).create_dict_fromDF(blacklist)

    # process each source, store the average array
    matrix_average = []
    for source in sources:
        # ready to extract region and write to output
        ofile = open(dirname+'/'+option.prefix+'_%s.csv' %(source), 'w')
        writer = csv.writer(ofile, lineterminator='\n', delimiter='\t')
        error = []
        isFirst = True

        # if need percentile coverage, load uniq value percentile here
        if option.normalize == 'percentile':
            ## load uniq values and counts
            genome_size = 0
            storeCounts = defaultdict(lambda: 0)
            for chr in nad.chrKeys[option.specie]:
                array = hdf5.load_chr(source, chr)
                genome_size += array.size
                uniq, counts = np.unique(array, return_counts=True)
                for index in range(len(uniq)):
                    storeCounts[uniq[index]] += counts[index]
            ## convert to percentile
            storeRank = {}
            startRank = 0
            endRank = 0
            for key, value in storeCounts.items():
                startRank = endRank + 1
                endRank += value
                rank = (startRank + endRank) / 2
                storeRank[key] = rank / genome_size
        elif option.normalize=='zscore':
            zmean = hdf.get_mean(source)
            zstd = np.sqrt(hdf.get_var(source))

        # extract coverage per region, store them, skip the region if overlap with blacklist
        matrix = []
        ids = []
        writer.writerow(["#Source: %s" %(source)])
        writer.writerow(["#Extract region: %s" %(option.region)])
        writer.writerow(["#Warning: There will be some regions disgarded because of region index out of bound"])
        writer.writerow(["#Error ID will be listed at the end of file"])
        size = region.loc[0, 'end'] - region.loc[0, 'start'] + 1
        for i in region.index:
            try:
                if (not blacklist is None and sum(blacklist.get_range(region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end']))==0) or (blacklist is None):
                    piece = hdf5.get_slice(source, region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end'])
                    if piece.size != size:
                        print("Skip region with length != the first region, ID: %s" %(region.loc[i, 'ID']))
                        continue
                    if region.loc[i,'strand'] == "-":
                        piece = piece[::-1]
                    ## convert to percentile if needed
                    if option.normalize=='percentile':
                        piece = toPercentile(piece, storeRank)
                    elif option.normalize=='zscore':
                        piece = (piece-zmean)/zstd
                    ofile.write("%s\t" %(region.loc[i, 'ID']))
                    writer.writerow(piece)
                    matrix.append(piece)
                    ids.append(region.loc[i, 'ID'])
            except Exception as e:
                error.append(region.loc[i, 'ID'])
        ofile.write("#Error id: ")
        writer.writerow(error)
        ofile.close()
        print("Error ID:\n%s" %(str(error)))

        # plot heatmap according to matrix
        sort_index = np.argsort(np.sum(matrix, axis=1))[::-1]
        matrix = np.vstack(matrix)[sort_index]   # sort by rowSum
        ids = np.array(ids)[sort_index]
        plt.figure(figsize=(12,12))
        plt.imshow(matrix, cmap='hot', interpolation='nearest', aspect='auto')
        plt.colorbar()
        plt.axis('off')
        plt.title("Coverage Heatmap\n%s in %s" %(source, option.region), fontsize=25, pad=10)
        plt.savefig(dirname+'/'+option.prefix+'_%s_heatmap.png' %(source), dpi=100)

        # save the id of each line in matrix
        with open(dirname + '/' + option.prefix + '_%s_id.txt' %(source), 'w') as f:
            for id in ids:
                f.write("%s\n" %id)

        # store the average of coverage
        array_average = np.average(matrix, axis=0)
        with open(dirname+'/'+option.prefix+'_%s_average.csv' %(source), 'w') as f:
            writer = csv.writer(f, lineterminator='\n', delimiter='\t')
            writer.writerow(["#Source: %s" %(source)])
            writer.writerow(["#Extract region: %s" %(option.region)])
            writer.writerow(["#Warning: There will be some regions disgarded because of region index out of bound"])
            writer.writerow(["#Error ID will be listed at the end of file"])
            writer.writerow(array_average)
            f.write("#Error id: ")
            writer.writerow(error)

        # plot composite plot
        ## First smooth the array
        array_average = np.convolve(array_average, np.ones((20,))/20, mode='valid')
        ## Then plot the composite plot
        plt.figure(figsize=(12,12))
        plt.plot(np.arange(array_average.size)-int(array_average.size/2), array_average, label=source)
        plt.legend(fontsize=15)
        plt.xlabel("Distance to Midpoint", fontsize=20)
        plt.ylabel("Average of Coverage Per Region", fontsize=20)
        plt.title("Composite Plot\n %s in %s" %(source, option.region), fontsize=25, pad=10)
        plt.tick_params(axis='both', labelsize=15)
        plt.savefig(dirname+'/'+option.prefix+'_%s_composite.png' %(source), dpi=100)
        ## Append to the average matrix
        matrix_average.append(array_average)
    # Return average matrix
    matrix_average = np.vstack(matrix_average).astype('float64')    
    return matrix_average

# load region in root processor, broadcast to workers
if rank == 0:
    region = pd.read_csv(option.region, header=None, sep='\t', comment='#', usecols=[0,1,2,3,4])
    region = region.astype({1:'int', 2:'int'})
    region.columns = ['chr', 'start', 'end', 'strand', 'ID']
else:
    region = None
region = comm.bcast(region, root=0)

# load blacklist if supplied, broadcast to workers
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
    interval = np.linspace(0, len(sources), size+1, dtype='int')
else:
    interval = None
interval = comm.bcast(interval, root=0)

# Time to work
## Open HDF5 object for each processor
hdf5 = nad.hdf5(option.hdf5, specie=option.specie, mpi=True, mpi_comm=MPI.COMM_WORLD)
## Process each source
matrix = get_coverage(hdf5, sources[interval[rank]:interval[rank+1]], region, blacklist)
collen = matrix.shape[1]
## Determine gather parameters
if rank == 0:
    recvbuf = np.empty([len(sources), collen], dtype='float64')
    split_sizes = []
    for  i in range(len(interval)-1):
        split_sizes = np.append(split_sizes, interval[i+1]-interval[i])
    split_sizes_output = split_sizes * collen
    displacement_output = np.insert(np.cumsum(split_sizes_output), 0, 0)[0:-1]
else:
    split_sizes_output = None
    displacement_output = None
    recvbuf = None
## Send average matrix to main processor
comm.Gatherv(sendbuf=matrix, recvbuf=[recvbuf, split_sizes_output, displacement_output, MPI.DOUBLE], root=0)
hdf5.close()

# Write the gathered matrix to output 
if rank == 0:
    df_average = pd.DataFrame(recvbuf, index=sources)
    ## save matrix
    df_average.to_csv(dirname+'/'+option.prefix+'_average.matrix', header=False, index=True, sep='\t')
    ## Plot composite plot
    plt.figure(figsize=(12,12))
    for index in range(len(sources)):
        plt.plot(np.arange(collen)-int(collen/2), recvbuf[index, :], label=sources[index])
    plt.legend(fontsize=15)
    plt.xlabel("Distance to Midpoint of Region", fontsize=20)
    plt.ylabel("Average of Coverage Per Region", fontsize=20)
    plt.title("Composite Plot in %s" %(option.region), fontsize=25, pad=10)
    plt.tick_params(axis='both', labelsize=15)
    plt.savefig(dirname+'/'+option.prefix+'_composite.png', dpi=100)








