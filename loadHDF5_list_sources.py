#!python3
#-*-coding:utf-8-*-

description = """loadHDF5_list_sources.py

List sources contained in HDF5 file

Example:
    $ python loadHDF5_list_sources.py mouse cache.h5

Dependencies:
    h5py

@Copyright 2019, Jianyu Yang, Southern Medical University
"""

import sys
import h5py

def main():
    filename = sys.argv[2]
    specie = sys.argv[1]

    hdf5 = h5py.File(filename, 'a')
    chrkey = list(hdf5['%s' %specie].keys())[0]
    print("Sources of %s in %s" %(specie, filename))
    for source in hdf5['%s/%s' %(specie, chrkey)].keys():
        print("\t%s" % source)

if __name__ == '__main__':
    main()