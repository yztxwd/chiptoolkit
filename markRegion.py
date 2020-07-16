#!python3

description = """

    Mark the region overlapped with the designated source in h5 database

    Example usage:
        python markRegion.py [specie] [database] [source] [region] [output]

    Requirements:
        first 3 columns of region file:
            chr start   end

"""

import sys, re
import numpy as np
import pandas as pd
import numpyArrayDict as nad

def main():
    specie, database, source, region, output = sys.argv[1:]

    # load the database
    h5 = nad.hdf5(database, specie=specie)

    # check each line of region
    with open(region, 'r') as f, open(output, 'w') as o:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                splits = line.strip().split("\t")

            chrom = "chr" + re.sub(r"chr", "", splits[0], flags=re.I)
            start = int(splits[1])
            end = int(splits[2])

            if sum(h5.get_slice(source, chrom, start, end)) > 0:
                newline = line.strip() + "\t1"
            else:
                newline = line.strip() + "\t0"

            o.write(newline + "\n")        

if __name__ == "__main__":
    if len(sys.argv) != 6 or sys.argv[1] in ['-h', '--help']:
        print(description)
        raise Exception("Exactly 4 parameters are required!")

    main()    
