#!python3.7

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from optparse import OptionParser, IndentedHelpFormatter

class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description

description = """

Plot heatmap and composite plot based on coverage matrix

Example usage:
    $ python heatmap.py -m matrix.csv -i -p matrix

Format Requirements:
    matrix:
        a tab delimited matrix describing the coverage per base per region

Dependencies:
    numpy
    pandas
    matplotlib

@Copyright 2020, Jianyu Yang, Southern Medical University

"""

parser = OptionParser(description= description, formatter = CustomHelpFormatter())
parser.add_option('-m', '--matrix', dest='matrix', help='Coverage matrix')
parser.add_option('-i', '--ignore', dest='ignore', action='store_true', help='Ignore the first column, set it as row label')
parser.add_option('-p', '--prefix', dest='prefix', help='Prefix of the output figure')
parser.add_option('-b', '--ylim_bottom_zero', dest='bottom', action='store_true', help='Set ylim bottom as zero')
parser.add_option('-s', '--scientific', dest='scientific', action='store_true', help='Use scientific notation')
option, argument = parser.parse_args()

# Load arguments
matrixFile = option.matrix
ignoreFirst = option.ignore
prefix = option.prefix
bottom = option.bottom
scientificNotation = option.scientific

# Load matrix info
matrixDF = pd.read_table(matrixFile, header=None, comment='#', index_col=0 if ignoreFirst else False)

# Generate matplotlib object
fig, axs = plt.subplots(2, 1, gridspec_kw={"height_ratios": [3,1]})

# Plot heatmap
## sort the matrix according to the sum of coverage per row
matrix = matrixDF.to_numpy()
matrixSort = matrix[np.argsort(matrix.sum(axis=1)),]
## log(n+1) transformation to ensure visibility of small value
matrixLog = np.log(matrixSort + 1)
## get heatmap
axs[0].imshow(matrixLog, cmap='hot', interpolation='nearest', aspect='auto')
axs[0].axis('off')

# Plot composite plot
## get the sum of coverage per base
matrixSum = np.sum(matrix, axis=0)
## smooth the array
matrixSmooth = np.convolve(matrixSum, np.ones((20,))/20, mode='valid')
## get composite plot
axs[1].plot(np.arange(matrixSmooth.size)-int(matrixSmooth.size/2), matrixSmooth)
axs[1].set_ylabel("Sum of coverage")
if bottom:
    axs[1].set_ylim(bottom=0)
if scientificNotation:
    axs[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# Save figure
fig.suptitle(prefix)
fig.savefig(prefix + ".png", dpi=100)

