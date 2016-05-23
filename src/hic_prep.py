import numpy as np
import pandas as pd
import sys
import argparse
import math

def lenstr2int(lenstr):
    # This function convers human friendly 
    # notation of length string into int value
    if(lenstr[-1:] in {'k', 'K'}):
        return(int(lenstr[:-1]) * 1000)
    elif(lenstr[-1:] in {'m', 'M'}):
        return(int(lenstr[:-1]) * 1000000)
    elif(lenstr[-1:] in {'g', 'G'}):
        return(int(lenstr[:-1]) * 1000000000)
    else:
        return(int(lenstr))

def bed_set(bedfile, chromosome, res):
    bed_set = set([])
    bed_all = pd.read_table(bedfile, header=None)
    bed_chr = bed_all[bed_all[0] == chromosome]    
    for interval in np.array(bed_chr[[1, 2]], dtype = np.int_):
        for i in range(int(interval[0] / res),
                       int(interval[1] / res) + 1):
            bed_set.add(i)
    return(bed_set)

def main_sub(infile, outfile, chromosome, 
             res, minsize, maxsize, 
             bedfile, normfile, expfile,
             log, distnorm):

    # read files

    if(bedfile != None):
        bed = bed_set(bedfile, chromosome, res)
    if(normfile != None):
        norm = np.loadtxt(normfile)
    if(expfile != None):
        exp  = np.loadtxt(expfile)

    data_obs = pd.read_table(infile, 
                            names = ['i', 'j', 'mij_raw'])

    # extract relevant data points

    i = np.array(data_obs['i'])
    j = np.array(data_obs['j'])
    dist = abs(i - j)

    if(bedfile != None):
        data_filtered = data_obs[(minsize <= dist) & 
                                 (dist <= maxsize) & 
                                 np.array([int(pos / res) in bed for pos in i]) &
                                 np.array([int(pos / res) in bed for pos in j])]
    else:
        data_filtered = data_obs[(minsize <= dist) & 
                                 (dist <= maxsize)]

    # convert raw observed contact frequencies

    if(normfile != None and expfile != None):
        if(log):
            log2_all = data_filtered.apply(lambda x: np.log(1.0 * x[2] / 
                                                            ((norm[int(x[0] / res)]) * 
                                                             (norm[int(x[1] / res)]) * 
                                                             exp[int(abs(x[1] - x[0]) / res)])) / np.log(2), 
                                           axis = 1)            
            log2 = log2_all[np.isfinite(log2_all)]
            if(distnorm):
                mij = (log2 - np.mean(log2)) / np.std(log2)
            else:
                mij = log2
            mij.name = 'mij'
            results = pd.concat([data_filtered[np.isfinite(log2_all)][['i', 'j']], mij], 
                                axis=1, join_axes=[mij.index])
        else:
            expoe_all = data_filtered.apply(lambda x: 1.0 * x[2] / 
                                            ((norm[int(x[0] / res)]) *
                                             (norm[int(x[1] / res)]) * 
                                             exp[int(abs(x[1] - x[0]) / res)]),
                                            axis = 1)            
            expoe = expoe_all[np.isfinite(expoe_all)]
            expoe.name = 'mij'
            results = pd.concat([data_filtered[np.isfinite(expoe_all)][['i', 'j']], expoe], 
                                axis=1, join_axes=[expoe.index])

    else:
        results = data

    if(outfile != None):
        results.to_csv(outfile, index = False, header = False, sep = '\t')


def main():
    parser = argparse.ArgumentParser(description='Hi-C prep')

    parser.add_argument('--version', action='version', version='%(prog)s 2016-05-23')

    parser.add_argument('-i', '--input',
                        metavar = 'i', default = None, required = True,
                        help = 'input file')

    parser.add_argument('-o', '--output',
                        metavar = 'o', default = None,
                        help = 'output file (default: stdout)')

    parser.add_argument('-c', '--chr',
                        metavar = 'c', default = None,
                        help = 'chromosome')

    parser.add_argument('-r', '--res',
                        metavar = 'r', default = None, required = True,
                        help = 'resolution')

    parser.add_argument('-m', '--min',
                        metavar = 'm', default = 0,
                        help = 'min (default = {})'.format(0))

    parser.add_argument('-M', '--max',
                        metavar = 'M', default = sys.maxsize,
                        help = 'Max (default = {})'.format(sys.maxsize))

    parser.add_argument('-b', '--bed',
                        metavar = 'b', default = None,
                        help = 'BED file to specify the region of interest')

    parser.add_argument('-n', '--norm',
                        metavar = 'n', default = None,
                        help = 'norm file')

    parser.add_argument('-e', '--exp',
                        metavar = 'e', default = None,
                        help = 'expected file')

    parser.add_argument('-l', '--log', action='store_true',
                        default = False,
                        help = 'log_2 conversion')

    parser.add_argument('-d', '--distnorm', action='store_true',
                        default = False,
                        help = 'convert distribution to N(0, 1)')


    args = parser.parse_args()

    

    print(args.input)
    print(args.output)
    print(args.chr)
    print(lenstr2int(args.res))
    print(lenstr2int(args.min))
    print(lenstr2int(args.max))
    print(args.bed)
    print(args.norm)
    print(args.exp)
    print(args.log)
    print(args.distnorm)


    main_sub(infile     = args.input, 
             outfile    = args.output,
             chromosome = args.chr, 
             res        = lenstr2int(args.res),
             minsize    = lenstr2int(args.min),
             maxsize    = lenstr2int(args.max), 
             bedfile    = args.bed,
             normfile   = args.norm,
             expfile    = args.exp,
             log        = args.log,
             distnorm   = args.distnorm)


if __name__ == '__main__':
    main()
