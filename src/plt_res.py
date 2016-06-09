import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import io
import math
import argparse

sns.set_style("ticks")
sns.set_context("paper", font_scale=2.0)

def read_res(file):
    names = ('axis', 'gamma', 'residuals', 'step_t', 'total_t')
    return(pd.read_table(file, names = names, skiprows = 1))

def plt_res(infile, outfile):
    data = read_res(infile)

    fig = plt.figure(figsize=(16,12))

    fig.suptitle(infile)

    ax1 = fig.add_subplot(2,2,1)
    ax1.scatter(range(len(data['gamma'])), data['gamma'])
    ax1.set_xlabel("iteration")
    ax1.set_xlim([0, len(data['gamma'])])
    ax1.set_ylim([1.2 * min(data['gamma']), 1.2 * max(data['gamma'])])
    ax1.set_ylabel("gamma")

    ax2 = fig.add_subplot(2,2,3)
    ax2.scatter(range(len(data['gamma'])), 
               data['residuals'])
    ax2.set_xlim([0, len(data['residuals'])])
    ax2.set_ylim([0.999 * min(data['residuals']), 1])
    ax2.set_xlabel("iteration")
    ax2.set_ylabel("residuals")

    ax3 = fig.add_subplot(2,2,4)
    ax3.scatter(range(len(data['gamma'])), 
                data['residuals'])
    ax3.set_xlim([0, len(data['residuals'])])
    ax3.set_ylim([0, 1])
    ax3.set_xlabel("iteration")
    ax3.set_ylabel("residuals")

    fig.savefig(outfile)

def main():
    parser = argparse.ArgumentParser(description='plot residuals curve')

    parser.add_argument('--version', action='version', version='%(prog)s 2016-06-03')
    parser.add_argument('-o', metavar = 'o', default = None, required = True,
                        help = 'output file')
    parser.add_argument('-i', metavar = 'i', default = None, required = True,
                        help = 'output file')
    args = parser.parse_args()

    plt_res(args.i, args.o)

if __name__ == '__main__':
    main()
