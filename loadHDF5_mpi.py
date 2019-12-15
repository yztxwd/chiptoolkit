import sys, os, csv
import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpi4py import MPI
from optparse import OptionParser, IndentedHelpFormatter
sys.path.append('/data2/yztxwd/scripts/python/Module')
import numpyArrayDict as nad


class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """loadHDF5_mpi.py

MPI version of loadHDF5.py
Load the coverage of sources in the given regions

Example:
    $ mpiexec -n 4 python loadHDF5.py --specie mouse --database cache.h5 --source zero,one
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
parser.add_option('-p', '--prefix',dest='prefix',default='loadHDF5',help='output file prefix')
parser.add_option('-o','--outDir',dest='outDir',default='.',help='output directory')
option, argument = parser.parse_args()

sources = option.source.split(",")
num = len(sources)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if num < size:
    raise Exception("Number of processors > Number of sources will cause unstability of program, please reduce core number")
    sys.exit(1)

# Global function definition
def get_coverage(hdf5, sources, region, blacklist):
    # get current working directory
    cwd = os.getcwd()+"/"

    # generate nad object if blacklist supplied
    if not blacklist is None:
        print("Use blacklist to exclude region")
        blacklist = nad.numpyArrayDict(specie=option.specie).create_dict_fromDF(blacklist)

    # make directory
    dirname = option.outDir
    if (not os.path.exists(dirname)):
        os.makedirs(dirname)

    for source in sources:
        # ready to extract region and write to output
        ofile = open(dirname+'/'+option.prefix+'_%s.csv' %(source), 'w')
        writer = csv.writer(ofile, lineterminator='\n', delimiter='\t')
        error = []
        isFirst = True

        # extract coverage per region, store them, skip the region if overlap with blacklist
        matrix = []
        writer.writerow(["#Source: %s" %(source)])
        writer.writerow(["#Extract region: %s" %(option.region)])
        writer.writerow(["#Warning: There will be some regions disgarded because of region index out of bound"])
        writer.writerow(["#Error ID will be listed at the end of file"])
        for i in region.index:
            try:
                if (not blacklist is None and sum(blacklist.get_range(region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end']))==0) or (blacklist is None):
                    piece = hdf5.get_slice(source, region.loc[i,'chr'], region.loc[i,'start'], region.loc[i,'end'])
                    if region.loc[i,'strand'] == "-":
                        piece = piece[::-1]
                    ofile.write("%s\t" %(region.loc[i, 'ID']))
                    writer.writerow(piece)
                    matrix.append(piece)
            except Exception as e:
                error.append(region.loc[i, 'ID'])
        ofile.write("#Error id: ")
        writer.writerow(error)
        ofile.close()
        print("Error ID:\n%s" %(str(error)))

        # plot heatmap according to matrix
        matrix = np.vstack(matrix)
        plt.figure(figsize=(12,12))
        plt.imshow(np.log(matrix+1), cmap='hot', interpolation='nearest', aspect='auto')
        print(np.log(matrix+1).shape)
        plt.colorbar()
        plt.axis('off')
        plt.title("Coverage Heatmap\n%s in %s" %(source, option.region), fontsize=25, pad=10)
        plt.savefig(dirname+'/'+option.prefix+'_%s_heatmap.png' %(source), dpi=100)

        # store the average of coverage
        matrix_average = np.average(matrix, axis=0)
        with open(dirname+'/'+option.prefix+'_%s_average.csv' %(source), 'w') as f:
            writer = csv.writer(f, lineterminator='\n', delimiter='\t')
            writer.writerow(["#Source: %s" %(source)])
            writer.writerow(["#Extract region: %s" %(option.region)])
            writer.writerow(["#Warning: There will be some regions disgarded because of region index out of bound"])
            writer.writerow(["#Error ID will be listed at the end of file"])
            writer.writerow(matrix_average)
            f.write("#Error id: ")
            writer.writerow(error)

        # plot composite plot
        ## First smooth the array
        matrix_average = np.convolve(matrix_average, np.ones((20,))/20, mode='valid')
        ## Then plot the composite plot
        plt.figure(figsize=(12,12))
        plt.plot(np.arange(matrix_average.size)-int(matrix_average.size/2), matrix_average, label=source)
        plt.legend(fontsize=15)
        plt.xlabel("Distance to Midpoint", fontsize=20)
        plt.ylabel("Average of Coverage Per Region", fontsize=20)
        plt.title("Composite Plot\n %s in %s" %(source, option.region), fontsize=25, pad=10)
        plt.tick_params(axis='both', labelsize=15)
        plt.savefig(dirname+'/'+option.prefix+'_%s_composite.png' %(source), dpi=100)

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
get_coverage(hdf5, sources[interval[rank]:interval[rank+1]], region, blacklist)

comm.Barrier()
hdf5.close()







