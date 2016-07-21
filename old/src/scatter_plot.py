import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import io
from scipy import stats
import math
import argparse

sns.set_style("ticks")
sns.set_context("paper", font_scale=2.0)

def read_results(file):
    col_names = ['i', 'j', 'obs', 'pred']
    res = pd.read_table(file, names = col_names, skiprows = 1)
    return(res)

def cor_dump(data):
    # Spearman's R
    print("Spearman's R")
    rho, rho_pval = stats.spearmanr(data['obs'], data['pred'])
    print("rho = {}\tpval = {}".format(rho, rho_pval))

    # Pearson's R
    print("Pearson's R")
    r, r_pval = stats.pearsonr(data['obs'], data['pred'])
    print("R   = {}\tpval = {}".format(r, r_pval))

def save_fig(data, png_file):
    grid = sns.JointGrid(x=data['obs'], y=data['pred'], space=0, size=6, ratio=50)
    grid.plot_joint(plt.scatter, color="g")
    grid.plot_marginals(sns.rugplot, height=1, color="g")

    grid.savefig(png_file)

def main_sub(res_file, png_file):
    data = read_results(res_file)
    print("results are now on memory")
    cor_dump(data)
    if(png_file != None):
        save_fig(data, png_file)

def main():
    parser = argparse.ArgumentParser(description='scatter plot')
    parser.add_argument('--version', action='version', version='%(prog)s 2016-06-01')
    parser.add_argument('-i', metavar = 'i', 
                        default = None, required = True,
                        help = 'input file (**.cmp.txt)')
    parser.add_argument('-o', metavar = 'o', 
                        default = None, 
                        help = 'output file (**.png)')
    args = parser.parse_args()

    main_sub(res_file = args.i, 
             png_file = args.o)

if __name__ == '__main__':
    main()
