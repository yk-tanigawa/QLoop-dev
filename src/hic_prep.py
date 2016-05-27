import numpy as np
import pandas as pd
from itertools import groupby
import sys, argparse, gzip
import logging
import math

# Fasta IO
def fasta_iter(fasta_name):
    '''
    given a fasta file. yield tuples of header, sequence
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    https://www.biostars.org/p/710/
    '''
    if((fasta_name[-3:] == '.gz') or 
       (fasta_name[-5:] == '.gzip')):
        with gzip.open(fasta_name, 'rb') as f:
            data = (x[1] for x in groupby(f, lambda line: line.decode('utf-8')[0] == ">"))
            for header in data:
                header = header.__next__().decode('utf-8')[1:].strip()
                seq = "".join(s.decode('utf-8').strip() for s in data.__next__())
                yield(header, seq)
    else:
        with open(fasta_name) as f:
            # ditch the boolean (x[0]) and just keep the header or sequence since
            # we know they alternate.
            data = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in data:
                # drop the ">"
                header = header.__next__()[1:].strip()
                # join all sequence lines to one.
                seq = "".join(s.strip() for s in data.__next__())
                yield(header, seq)

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

def fasta_chk_gap(seq, res, margin):
    has_not_gap = np.array(['N' not in seq[i * res - margin : (i + 1) * res + margin].upper() for i in range(1 + int(len(seq) / res))])
    has_not_gap[0] = has_not_gap[-1] = False                                               
    return(has_not_gap)

def hic_prep(logger, 
             infile, outfile, chromosome, 
             res, margin, minsize, maxsize, 
             fastafile, bedfile, normfile, expfile,             
             log, distnorm):

    # read files
    if(fastafile != None):
        seqs = {}
        for (head, seq) in fasta_iter(fastafile):
            seqs[head] = seq
        fasta = fasta_chk_gap(seqs[chromosome], res, margin)
        logger.debug('Fasta file has been loaded')

    if(bedfile != None):
        bed = bed_set(bedfile, chromosome, res)
        logger.debug('BED file has been loaded')
    if(normfile != None):
        norm = np.loadtxt(normfile)
        logger.debug('normalization vector file has been loaded')
    if(expfile != None):
        exp  = np.loadtxt(expfile)
        logger.debug('expected value vector file has been loaded')

    data_obs = pd.read_table(infile, 
                            names = ['i', 'j', 'mij_raw'])
    logger.debug('raw data file has been loaded')

    # extract relevant data points

    i = np.array(data_obs['i'])
    j = np.array(data_obs['j'])
    dist = abs(i - j)

    if(fastafile != None):
        if(bedfile != None):
            data_filtered = data_obs[(minsize <= dist) & 
                                     (dist <= maxsize) & 
                                     np.array([fasta[int(pos / res)] for pos in i]) &
                                     np.array([fasta[int(pos / res)] for pos in j]) &
                                     np.array([int(pos / res) in bed for pos in i]) &
                                     np.array([int(pos / res) in bed for pos in j])]
        else:
            data_filtered = data_obs[(minsize <= dist) & 
                                     (dist <= maxsize) &
                                     np.array([fasta[int(pos / res)] for pos in i]) &
                                     np.array([fasta[int(pos / res)] for pos in j])]
    else:
        if(bedfile != None):
            data_filtered = data_obs[(minsize <= dist) & 
                                     (dist <= maxsize) & 
                                     np.array([int(pos / res) in bed for pos in i]) &
                                     np.array([int(pos / res) in bed for pos in j])]
        else:
            data_filtered = data_obs[(minsize <= dist) & 
                                     (dist <= maxsize)]
    logger.debug('relevant data points are extracted')

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

    logger.debug('data conversion is finished')

    if(outfile != None):
        logger.debug('writing the results to {}'.format(outfile))
        results.to_csv(outfile, index = False, header = False, sep = '\t')

def hic_prep_args_dump(logger, 
                       infile, outfile, chromosome, 
                       res, margin, minsize, maxsize, 
                       fastafile, bedfile, normfile, expfile,             
                       log, distnorm):
    logger.info('infile     : {}'.format(infile))
    logger.info('outfile    : {}'.format(outfile))
    logger.info('chromosome : {}'.format(chromosome))
    logger.info('resolution : {}'.format(res))
    logger.info('margin     : {}'.format(margin))
    logger.info('min size   : {}'.format(minsize))
    logger.info('max size   : {}'.format(maxsize))
    logger.info('fasta file : {}'.format(fastafile))
    logger.info('bed file   : {}'.format(bedfile))
    logger.info('norm. file : {}'.format(normfile))
    logger.info('exp. file  : {}'.format(expfile))
    logger.info('log conv.  : {}'.format(log))
    logger.info('normalize  : {}'.format(distnorm))
    

def hic_prep_main():
    ## argparse
    parser = argparse.ArgumentParser(description='Hi-C prep')

    parser.add_argument('--version', action='version', version='%(prog)s 2016-05-24')

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

    parser.add_argument('-a', '--margin',
                        metavar = 'a', default = None, 
                        help = 'margin')

    parser.add_argument('-m', '--min',
                        metavar = 'm', default = 0,
                        help = 'min (default = {})'.format(0))

    parser.add_argument('-M', '--max',
                        metavar = 'M', default = sys.maxsize,
                        help = 'Max (default = {})'.format(sys.maxsize))

    parser.add_argument('-f', '--fasta',
                        metavar = 'f', default = None,
                        help = 'Fasta file (to check gap)')

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

    parser.add_argument('-v', '--verbose', type = int,
                        metavar = 'v', default = logging.DEBUG,
                        help = 'verbose level (default : {})'.format(logging.DEBUG))

    args = parser.parse_args()

    ## logging
    # create logger
    logger = logging.getLogger('logger')
    logger.setLevel(args.verbose)

    # create console handler and set level to debug
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)

    # create formatter
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')

    # add formatter to ch
    ch.setFormatter(formatter)

    # add ch to logger
    logger.addHandler(ch)


    logger.debug('Hello. Starting the program.')
    
    hic_prep_args_dump(logger     = logger,
                       infile     = args.input, 
                       outfile    = args.output,
                       chromosome = args.chr, 
                       res        = lenstr2int(args.res),
                       margin     = lenstr2int(args.margin),
                       minsize    = lenstr2int(args.min),
                       maxsize    = lenstr2int(args.max), 
                       fastafile  = args.fasta,
                       bedfile    = args.bed,
                       normfile   = args.norm,
                       expfile    = args.exp,
                       log        = args.log,
                       distnorm   = args.distnorm)

    hic_prep(logger     = logger, 
             infile     = args.input, 
             outfile    = args.output,
             chromosome = args.chr, 
             res        = lenstr2int(args.res),
             margin     = lenstr2int(args.margin),
             minsize    = lenstr2int(args.min),
             maxsize    = lenstr2int(args.max), 
             fastafile  = args.fasta,
             bedfile    = args.bed,
             normfile   = args.norm,
             expfile    = args.exp,
             log        = args.log,
             distnorm   = args.distnorm)


if __name__ == '__main__':
    hic_prep_main()
