import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import io
import math
import argparse
import os.path

sns.set_style("ticks")
sns.set_context("paper", font_scale=2.0)

def read_res(file):
    names = ('axis', 'gamma', 'residuals', 'step_t', 'total_t')
    return(pd.read_table(file, names = names, skiprows = 1))

def set_plt_style(n):
    colors = sns.husl_palette(n)
    markers = ['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd'] * (int(n / 13) + 1)
    linestyles = ['-', ':', '-.', '--', ':', '-.', '--'] * (int(n / 7) + 1)
    return({'color': colors,
            'marker': markers,
            'linestyle': linestyles})

def plt_res(infiles, outfile, title):
    fig = plt.figure(figsize=(16,12))
    fig.suptitle(title, fontsize = 30)
    plt_res_sub(fig, infiles)
    fig.show()
    if(outfile != None):
        fig.savefig(outfile)

def plt_res_sub(fig, infiles):    
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)

    data = [read_res(infile) for infile in infiles]
    plt_style = set_plt_style(len(data))
    
    ckps_max = 0
    for i in range(len(data)):
    
        ax1.scatter(range(len(data[i]['gamma'])), data[i]['gamma'],
                    color = plt_style['color'][i],
                    marker = plt_style['marker'][i],
                    linestyle = plt_style['linestyle'][i])

        weight = data[i].groupby('axis').sum()['gamma']
        ax2.scatter(weight.index, weight,
                    color = plt_style['color'][i],
                    marker = plt_style['marker'][i],
                    linestyle = plt_style['linestyle'][i])
        if(ckps_max < max(weight.index)):
            ckps_max = max(weight.index)

        ax3.scatter(range(len(data[i]['gamma'])), data[i]['residuals'],
                    color = plt_style['color'][i],
                    marker = plt_style['marker'][i],
                    linestyle = plt_style['linestyle'][i])
        ax4.scatter(range(len(data[i]['gamma'])), data[i]['residuals'],
                    color = plt_style['color'][i],
                    marker = plt_style['marker'][i],
                    linestyle = plt_style['linestyle'][i],
                    label=os.path.basename(infiles[i]))
                    
    iter_num_max = max([len(data[i]['gamma']) for i in range(len(data))])
    ax1.set_xlim([0, iter_num_max])
    ax2.set_xlim([0, ckps_max])
    ax3.set_xlim([0, iter_num_max])
    ax4.set_xlim([0, iter_num_max])
    ax1.set_xlabel("iteration")
    ax2.set_xlabel("$k$-mer pair (id)")
    ax3.set_xlabel("iteration")
    ax4.set_xlabel("iteration")

    ax1.set_title("$\\nu\gamma$ in each iteration ")
    ax2.set_title("weight $w$ of $k$-mer pairs")
    ax3.set_title("residual curve")
    ax4.set_title("residual curve")

    ax1.yaxis.offsetText.set_fontsize(30)
    ax2.yaxis.offsetText.set_fontsize(30)
    ax1.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
    ax2.ticklabel_format(style='sci',axis='y',scilimits=(0,0))

    ax1.set_ylabel("$\\nu\gamma$", fontsize=30)                
    ax1.set_ylim([1.2 * min([min(data[i]['gamma']) for i in range(len(data))]),
                  1.2 * max([max(data[i]['gamma']) for i in range(len(data))])])

    ax2.set_ylabel("weight $w$")

    ax3.set_ylabel("residuals")    
    ax3.set_ylim([0.99 * min([min(data[i]['residuals']) for i in range(len(data))]), 1])

    ax4.set_ylabel("residuals")
    ax4.set_ylim([0, 1])

    ax4.legend(loc='lower left', ncol=1)

def plt_res_old(infile, outfile):
    data = read_res(infile)

    fig = plt.figure(figsize=(16,12))

    fig.suptitle(infile)

    ax1 = fig.add_subplot(2,2,1)
    ax1.scatter(range(len(data['gamma'])), data['gamma'])
    ax1.set_xlabel("iteration")
    ax1.set_xlim([0, len(data['gamma'])])
    ax1.set_ylim([1.2 * min(data['gamma']), 1.2 * max(data['gamma'])])
    ax1.set_ylabel("gamma")

    ax3 = fig.add_subplot(2,2,3)
    ax3.scatter(range(len(data['gamma'])), 
               data['residuals'])
    ax3.set_xlim([0, len(data['residuals'])])
    ax3.set_ylim([0.999 * min(data['residuals']), 1])
    ax3.set_xlabel("iteration")
    ax3.set_ylabel("residuals")

    ax4 = fig.add_subplot(2,2,4)
    ax4.scatter(range(len(data['gamma'])), 
                data['residuals'])
    ax4.set_xlim([0, len(data['residuals'])])
    ax4.set_ylim([0, 1])
    ax4.set_xlabel("iteration")
    ax4.set_ylabel("residuals")

    fig.savefig(outfile)

def main():
    parser = argparse.ArgumentParser(description='plot residuals curve')

    parser.add_argument('--version', action='version', version='%(prog)s 2016-06-12')
    parser.add_argument('-o', metavar = 'o', default = None, 
                        help = 'output file')
    parser.add_argument('-i', metavar = 'i', default = None, required = True,
                        nargs='+',
                        help = 'input file')
    parser.add_argument('-t', metavar = 't', default = "summary of the results", 
                        help = 'figure title')
    args = parser.parse_args()

    plt_res(infiles = args.i, outfile = args.o, title = args.t)

if __name__ == '__main__':
    main()
