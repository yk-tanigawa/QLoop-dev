import seaborn as sns
import numpy as np
import pandas as pd
from scipy import io
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
    return(pd.read_table(file, names = names, skiprows = 1))


def res_to_merged_sub(res_file, ckps_file):
    '''
    This function reads both canonical k-mer file results file and merge them into one
    '''
    (ckps_id, ckps) = read_ckps(ckps_file)
    res = read_res(res_file)
    merged = pd.concat([ckps.loc[res.axis],  ckps_id.loc[res.axis]], axis = 1)
    merged = merged.reset_index(drop = True)
    merged = pd.concat([merged, res], axis = 1)
    return(merged.loc[1:])

def res_to_merged(res_file, ckps_file, out_file):
    if(out_file == None):
        out_file = res_file + '.merged'

    merged = res_to_merged_sub(res_file, ckps_file)
    merged.to_csv(out_file)


def main():
    parser = argparse.ArgumentParser(description='a converter from results file to merged results')

    parser.add_argument('--version', action='version', 
                        version='%(prog)s 2016-06-08')
    parser.add_argument('-r', metavar = 'r', 
                        default = None, required = True,
                        help = 'results file')
    parser.add_argument('-c', metavar = 'c', 
                        default = None, required = True,
                        help = 'canonical k-mer pair list file')
    parser.add_argument('-o', metavar = 'o', 
                        default = None, 
                        help = 'output file')
    args = parser.parse_args()

    res_to_merged(res_file  = args.r,
                  ckps_file = args.c, 
                  out_file  = args.o)

if __name__ == '__main__':
    main()
