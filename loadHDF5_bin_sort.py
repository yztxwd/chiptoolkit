#!/usr/bin/python3.6
#coding=utf-8

"""
This script aims to chop the genome into small bins, then sort them by reads number
"""

import pandas as pd
import numpy as np
import sys, csv
sys.path.append("/data2/yztxwd/scripts/python/Module")
import numpyArrayDict as nad

def chop_bin_sum_norm(database, specie, source, bin_size = 1000):
    # open hdf5 file
    hdf5 = nad.hdf5(database, specie = specie)
    
    # load numpyArrayDict
    mid = hdf5.load_numpyArrayDict(source).store
    
    # construct a dictionary containing bins per chromosome
    store = {}
    for i in mid.keys():
        extra = mid[i].size%bin_size
        if(extra!=0):
            append = bin_size - extra
            mid[i] = np.append(mid[i], np.zeros(append))
        store[i] = mid[i].reshape([-1, bin_size]).sum(axis=1)
        print(store[i].shape)
    # concatenate all bins
    all_chr = np.concatenate(list(store.values()))
    print(all_chr.shape)
    all_chr = all_chr[all_chr.argsort()]
    all_chr_norm = all_chr / all_chr.sum()
    cum_sum = all_chr_norm.cumsum()

    return all_chr, all_chr_norm, cum_sum
    
if __name__ == "__main__":
    if(sys.argv[1] in ["-h", "--help"]):
        print("""
        This script aims to chop the genome into small bins, then sort them by normalized reads number

        example usage:
            python script.py hdf5_filename specie source bin_size prefix ifZero
            python script.py ~/Data/CACHE.hdf5 yeast index 1000 test
        """)
        sys.exit(1)
    db, specie, source, bin_size, prefix = sys.argv[1:6]
    bin_size = int(bin_size)
    bin_sum, bin_norm, bin_norm_cum_sum = chop_bin_sum_norm(db, specie, source, bin_size)
    ofile = open(prefix+"_bin_sum.txt","w")
    writer = csv.writer(ofile, lineterminator="\n", delimiter="\t")
    
    # write to output
    writer.writerow(["#Source: %s" %(source)])
    writer.writerow(["#Bin size: %s" %(bin_size)])
    writer.writerow(["bin_counts","normalized_bin_counts","normalized_bin_counts_cumsum"])
    for i in range(0, bin_norm.size):
        writer.writerow([bin_sum[i], bin_norm[i], bin_norm_cum_sum[i]])
    ofile.close()


    
