#!python3.7

import numpy as np
import pandas as pd

from optparse import OptionParser, IndentedHelpFormatter

class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """

Plot heatmap and composite plot based on coverage matrix

Example


"""

parser = OptionParser(description= description, formatter = CustomHelpFormatter())
parser.add_option('-m', '--matrix', dest='matrix', help='Coverage matrix')
parser.add_option('-i', '--ignore', dest='ignore', action='store_true', help='Ignore the first column, set it as row label')
parser.add_option('-p', '--prefix', dest='prefix', help='Prefix of the output figure')
option, argument = parser.parse_args()

# Load arguments
matrixFile = option.matrix
ignoreFirst = option.ignore

# Load matrix info
matrixDF = pd.read_table(matrixFile, header=None, comment='#', index_col=0 if ignoreFirst else False)



# Plot heatmap
## sort the matrix according to the sum of coverage per row
matrix = matrixDF.to_numpy()
matrixSort = matrix[np.argsort(matrix.sum(axis=1)),]
## log(n+1) transformation to ensure visibility of small value
matrixLog = np.log(matrixSort + 1)
## get heatmap

# Plot composite plot
## get the sum of coverage per base
matrixSum = np.sum(matrix, axis=0)
## smooth the array
matrixSmooth = np.convolve(matrixSum, np.ones((20,))/20, mode='valid')
## get composite plot

