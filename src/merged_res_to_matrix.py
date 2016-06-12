import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import io
from scipy.sparse import dok_matrix
import math
import argparse

sns.set_style("ticks")
sns.set_context("paper", font_scale=2.0)

def read_merged_results(file):
    '''
    read merged results file (tsv format)
    '''
    names = ['iter',
             'kmer1','kmer2','revcmp1','revcmp2',
             'kmer1_id','kmer2_id','revcmp1_id','revcmp2_id',
             'axis','gamma','residuals',
             'step_t','total_t']
    return(pd.read_table(file, 
                         names = names, 
                         sep = ',', 
                         skiprows = 1, 
                         index_col = names[0]))

def kmer_interaction_mat(data, col = 'gamma'):
    kmers = sorted(list(set([kmer for kmer in data.kmer1] + 
                            [kmer for kmer in data.kmer2])))
    kmer_ids = dict(zip(kmers, range(len(kmers))))
    
    mat = dok_matrix((len(kmers), len(kmers)), dtype = np.int_)
    for i in data.index:
        mat[kmer_ids[data.loc[i, 'kmer1']], 
            kmer_ids[data.loc[i, 'kmer2']]] += data.loc[i, col]
    return(mat, kmers)

def merged_res_to_matrix(infile, outfile, col = 'output'):
    if(outfile == None):
        outfile = infile + '.top5' + '.csv'
    
    merged_res = read_merged_results(infile)
    
    top5 = merged_res['gamma'] >= sorted(merged_res['gamma'])[int(len(merged_res) * 0.95)]

    merged_res[col] = round(100 * merged_res['gamma'] / min(merged_res[top5]['gamma']))

    (mat, kmers) = kmer_interaction_mat(merged_res[top5], col = col)

    dense_mat = pd.DataFrame(data = mat.todense(), 
                             index = kmers, columns = kmers)

    dense_mat.to_csv(outfile, sep = ' ')

def main():
    parser = argparse.ArgumentParser(description='a converter from marged results to csv file (to plot circos tableviewer plot)')

    parser.add_argument('--version', action='version', 
                        version='%(prog)s 2016-06-08')
    parser.add_argument('-i', metavar = 'i', 
                        default = None, required = True,
                        help = 'input file')
    parser.add_argument('-o', metavar = 'o', 
                        default = None,
                        help = 'output file')
    args = parser.parse_args()

    merged_res_to_matrix(infile = args.i, outfile = args.o)

if __name__ == '__main__':
    main()
