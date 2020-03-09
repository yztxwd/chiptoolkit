#!/usr/bin/python python3.6
#coding = utf-8

import re
import pandas as pd
import numpy as np
import h5py

description="""
    This module is used to store genome wide information on each base pair (i.e., depth), it employs a 0-based coordinate system (which is the default of numpy)
    When dump into or load from this module, you should use coordinate under 1-based coordinate system to ensure it works in the right way
"""

chrKeys = {
    "mouse": [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 
        'chrX', 'chrY', 'chrM'
    ],
    "yeast": [
        "chrI", "chrII", "chrIII", "chrIV", "chrIX", "chrM", "chrV", "chrVI", "chrVII",
        "chrVIII", "chrX", "chrXI", "chrXII", "chrXIII", "chrXIV", "chrXV", "chrXVI"
    ]
}

chrSize = {
    "mouse": [
        195471971, 182113224, 160039680, 156508116, 151834684, 149736546, 145441459, 129401213, 
        124595110, 130694993, 122082543, 120129022, 120421639, 124902244, 104043685, 98207768, 
        94987271, 90702639, 61431566, 171031299, 91744698, 16299
    ],
    "yeast": [
        230218, 813184, 316620, 1531933, 439888, 85779, 576874, 270161, 1090940, 562643, 
        745751, 666816, 1078177, 924431, 784333, 1091291, 948066
    ]
}

class numpyArrayDict(object):
    """
    Standard numpyArrayDic object, 1-based
    format:
    {
    'chr1': np.array([0, 1, 2, 3, 5, 66, 11])
    'chr2': np.array([0, 1, 2, 3, 5, 66, 11])
    ...
    'chrM': np.array([3, 5, 6, 7, 8, 11, 23])
    }
    """
    def __init__(self, specie='mouse'):
        self.store = {}
        self.specie = specie
        self.chrKeys = chrKeys
        self.chrSize = chrSize

    def normalize_dict(self, depth):
        """
        normalize self.store to 10M reads, each row counts as one read
        :return: normalized store
        """
        for i in self.store.keys():
            self.store[i] *= 10**7/depth

    def generate_dict(self):
        """
        Generate a dictionary containing coverage on mouse genome, [datatype: float16]
        output dictionary:
        key1: numpy.array[1,23,4,5,6,...]
        key2: numpy.array[1,26,7,8,9,...]
        """

        # Generate a dictionary containing the per base depth of all chromosomes
        for i in self.chrKeys[self.specie]:
            self.store[i] = np.zeros(self.chrSize[self.specie][self.chrKeys[self.specie].index(i)], dtype='float32')

    def fill_dict(self, chunk):
        """
        Fill in values in store according to the dataframe chunk
        """
        # Standardize the chromosome name to avoid ambiguous names of the same chromosome, i.e., chr1, Chr1, 1...
        chunk['chr'] = chunk.apply(lambda row: "chr" + re.sub(r"chr", "", row['chr'], flags=re.I), axis=1)
        
        def func(chr, start, end, depth):
            try:
                self.store[chr][(start-1):(end)] += depth
                return 0
            except KeyError:
                print("KeyError: %s" %(chr))
                return 0
        vfunc = np.vectorize(func, otypes=[int])
        vfunc(chunk['chr'].values, chunk['start'].values, chunk['end'].values, chunk['depth'].values)

    def create_dict_fromfile(self, filename, chunksize=10**6):
        """
        Create dictionary according to the input file
        :param filename:
        format: chr start end depth (separated by '\t')
        :return: self
        """
        self.generate_dict()
        chunks = pd.read_csv(filename, header=None, sep='\t', chunksize=chunksize, comment="#", usecols=[0,1,2,3], 
                names=['chr', 'start', 'end', 'depth'], dtype={'chr':str, 'start':int, 'end':int, 'depth':float})
        for chunk in chunks:
            self.fill_dict(chunk)
        return self

    def create_dict_fromDF(self, dataframe):
        """
        Create dictionary according to the input pandas dataframe
        :param dataframe:
        format: chr start end depth
        :return: self
        """
        # format dataframe
        dataframe.columns = ['chr', 'start', 'end', 'depth']
        dataframe.astype({'chr':str, 'start':int, 'end':int, 'depth':float})

        # dump into numpy array
        self.generate_dict()
        self.fill_dict(dataframe)
        return self

    def get_dict(self):
        """
        Return whole dictionary
        """
        return self.store

    def get_range(self, chr, start, end):
        """
        Return piece of genome according to given start, end coordinate
        :param chr:
        :param start:
        :param end:
        :return: numpy array
        """
        start = start-1
        end = end
        if (start<0 or end>self.chrSize[self.specie][chrKeys[self.specie].index(chr)]):
            raise Exception('Index out of chromosome bound!')
        return self.store[chr][start:end]

    def get_point(self, chr, index):
        index = index-1
        if (index<0 or index>=self.chrSize[self.specie][chrKeys[self.specie].index(chr)]):
            raise Exception('Index out of chromosome bound!')
        return self.store[chr][index]

    def get_size(self, chr):
        return self.chrSize[self.specie][self.chrKeys[self.specie].index(chr)]

    def dump_hdf5(self, hdf5, sourceName):
        """
        Save numpy dictionary to hdf5 file
        :param hdf5:
        :return:
        """
        file = h5py.File(hdf5)

        # check if all keys exist in hdf file, create if not
        if not (self.specie in file.keys()):
            file.create_group(self.specie)
        for i in self.chrKeys[self.specie]:
            if not (i in file[self.specie].keys()):
                file[self.specie].create_group(i)

        # dump numpy array dictionary into hdf5
        for i in self.chrKeys[self.specie]:
            file['/%s/%s' %(self.specie, i)].create_dataset(sourceName, data=self.store[i], compression='gzip')

    # deprecated
    def add_dict_fromfile(self, filename, addHeader=True, chunksize=10**6):
        chunks = pd.read_csv(filename, header=None, sep='\t', chunksize=chunksize, comment="#", usecols=[0,1,2,3])
        for chunk in chunks:
            self.fill_dict(chunk, addHeader=addHeader)

class hdf5(object):
    """
    handle hdf5 file containing numpy array dictionary
    """
    def __init__(self, hdf5, specie='mouse'):
        self.hdf5 = h5py.File(hdf5)
        self.specie = specie
        self.chrKeys = chrKeys
        self.chrSize = chrSize

    def get_point(self, sourceName, chr, coor):
        zero_based_coor = coor - 1
        if (zero_based_coor<0 or zero_based_coor>=self.chrSize[self.specie][chrKeys[self.specie].index(chr)]):
            raise Exception('Index out of chromosome bound!')
        returnValue = self.hdf5['/%s/%s/%s' %(self.specie, chr, sourceName)][zero_based_coor]
        return returnValue

    def get_slice(self, sourceName, chr, start, end):
        start = start-1
        end = end
        if (start<0 or end>self.chrSize[self.specie][chrKeys[self.specie].index(chr)]):
            raise Exception('Index out of chromosome bound!')
        returnArray = self.hdf5['/%s/%s/%s' %(self.specie, chr, sourceName)][start:end]
        return returnArray

    def load_numpyArrayDict(self, sourceName):
        """
        Load entire numpyArrayDict class from hdf5 file
        :param sourceName:
        :return:
        """
        store = {}
        for i in self.chrKeys[self.specie]:
            store[i] = self.hdf5['/%s/%s/%s' %(self.specie, i, sourceName)][:]

        returnClass = numpyArrayDict(specie=self.specie)
        returnClass.store = store
        return returnClass

    def list_source(self):
        """
        List names of all source in hdf5 file
        :return:
        """
        source_list = []
        for i in self.hdf5['/%s/%s' %(self.specie, self.chrKeys[self.specie][0])].keys():
            source_list.append(i)
        return source_list

    def delete_source(self, sourceName):
        """
        Delete dataset with sourceName
        :return:
        """
        try:
            for i in self.chrKeys[self.specie]:
                del self.hdf5['/%s/%s/%s' %(self.specie, i, sourceName)]
        except KeyError:
            print("%s doesn't exist in hdf5 file!" %(sourceName))
    
    def close(self):
        self.hdf5.close()

def main():
    pass


def find_nearest(chr, coor, hdf5, source,specie="yeast", strand="+"):
    range = 0
    sum = 0

    while(sum==0):
        range += 1000
        start = coor-range
        end = coor+range
        try:
            reg = hdf5.get_slice(source, chr, start, end)  #Extract +/- 1000 bp region around dyad to get nearest reference dyad locations
        except:
            if(start<=0):
                start = 1
                reg = hdf5.get_slice(source, chr, start, end)
                append_size = 2*range + 1 - reg.size
                reg = np.append(np.zeros(append_size), reg)
            else:
                end = numpyArrayDict(specie=specie).get_size(chr)
                reg = hdf5.get_slice(source, chr, start, end)
                append_size = 2*range + 1 - reg.size
                reg = np.append(reg, np.zeros(append_size))
            
        sum = reg.sum()

    #Get the relative positions of reference dyads
    index_dyad = np.where(reg > 0)[0]  
    #Compute distance to each reference dyads in region
    distance = index_dyad - (reg.size - range - 1)
    #Get the nearest reference dyad
    if (strand=="+"):
        distance = distance[distance>0]
    elif (strand=="-"):
        distance = distance[distance<0]
    min_dist_index = np.argmin(np.abs(distance))
    min_dist_dyad = coor + distance[min_dist_index]

    return min_dist_dyad


if __name__ == '__main__':
    main()
