# chiptoolkit

A collection of scripts for retrieving Genomics coverage and downstream analysis, this project is based on HDF5 to improve performance when working on large numbers of bigwig genomic tracks.

- **dumpHDF5.py**

    Compute per base coverage from a 1-based bedgraph format file, then dump into the specified hdf5 file

- **loadHDF5.py**

    Load coverage in the regions provided, output the average coverage per base, heatmap and composite plot. **All regions must have the same length**

- **loadHDF5_mpi.py**

    MPI version of loadHDF5.py

- **loadHDF5_pearson.py**

    Load genome-wide coverage or coverage in the regions if provided region file, concatenate them and chop into bins, compute pearson/spearman correlation of these bins, output correlation matrix, heatmap.

- **loadHDF5_pearson_mpi.py**

    MPI version of loadHDF5_pearson_mpi.py
