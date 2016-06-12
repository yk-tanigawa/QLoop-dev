import numpy as np
import pandas as pd
from scipy import io
from scipy.sparse import dok_matrix
import os.path
import math, argparse

def read_ckps(file):
    '''
    This function reads a canonical k-mer file
    '''
    head_id = ['kmer1_id', 'kmer2_id', 'revcmp1_id', 'revcmp2_id']
    header  = ['kmer1', 'kmer2', 'revcmp1', 'revcmp2']
    ckps_all = pd.read_table(file, names = head_id + header)
    ckps_id = ckps_all[head_id]
    ckps = ckps_all[header]    
    return(ckps_id, ckps)

def read_res(file):
    '''
    This function reads results file from L2 Boosting
    '''
    names = ('axis', 'gamma', 'residuals', 'step_t', 'total_t')
    return(pd.read_table(file, names = names, skiprows = 1)[1:])

def top(percentile, df):
    return(df >= sorted(df)[int(len(df) * percentile)])

def kmer_interaction_mat(data, col = 'gamma'):
    '''
    This function converts a data frame into a sparse matrix
    (a preparation step for circos plot)
    '''
    kmers = sorted(list(set([kmer for kmer in data.kmer1] + 
                            [kmer for kmer in data.kmer2])))
    kmer_ids = dict(zip(kmers, range(len(kmers))))
    
    mat = dok_matrix((len(kmers), len(kmers)), dtype = np.int_)
    for i in data.index:
        mat[kmer_ids[data.loc[i, 'kmer1']], 
            kmer_ids[data.loc[i, 'kmer2']]] += data.loc[i, col]
    return(mat, kmers)

def res_to_merged(res_file, ckps_file, out_file, percentile):
    # read a canonical k-mer pair list from a file
    (ckps_id, ckps) = read_ckps(ckps_file)
    # read the weights of k-mer pairs from the result file
    weight = read_res(res_file).groupby('axis').sum()['gamma']
    # create a unified table
    merged = ckps.loc[weight.index]
    merged['w'] = weight

    if(out_file == None):
        print(merged[top(percentile, merged['w'])].sort_values(by = 'w', ascending = False))
    else:
        # compute values for circos plot in an ad-hock manner
        col = 'circos'
        merged[col] = round(100 * merged['w'] / min(merged[top(percentile, merged['w'])]['w']))

        # convert selected results into a sparse matrix
        (mat, kmers) = kmer_interaction_mat(merged[top(percentile, merged['w'])], col = col)

        # convert the sparse matrix into a dense matrix (to generate csv file)
        dense_mat = pd.DataFrame(data = mat.todense(), 
                                 index = kmers, columns = kmers)

        # save the results as csv file
        dense_mat.to_csv(out_file, sep = ' ')

def main():
    parser = argparse.ArgumentParser(description='a converter from results file to merged results')

    parser.add_argument('--version', action='version', 
                        version='%(prog)s 2016-06-12')
    parser.add_argument('-r', metavar = 'r', 
                        default = None, required = True,
                        help = 'results file')
    parser.add_argument('-c', metavar = 'c', 
                        default = None, required = True,
                        help = 'canonical k-mer pair list file')
    parser.add_argument('-o', metavar = 'o', 
                        default = None, 
                        help = 'output file (csv)')
    parser.add_argument('-p', metavar = 'p', type=float,
                        default = 0.0, 
                        help = 'percentile')
    args = parser.parse_args()

    res_to_merged(res_file  = args.r,
                  ckps_file = args.c, 
                  out_file  = args.o,
                  percentile = args.p)

if __name__ == '__main__':
    main()
