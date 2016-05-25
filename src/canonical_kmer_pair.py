import argparse, itertools

def revcomp(mer, comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}):
    return("".join([comp[mer[i-1:i]] for i in range(len(mer), 0, -1)]))

def mer2bin(mer, table = {'A':'00', 'C':'01', 'G':'10', 'T':'11'}):
    return(int("".join([table[mer[i:i+1]] for i in range(len(mer))]), 2))

def canonical_kmer_pairs(k, eliminate = None, c = ['A', 'C', 'G', 'T']):
    if(eliminate != None):
        mers = ["".join(element) for element in itertools.product(c, repeat=k) if eliminate not in "".join(element)]
    else:
        mers = ["".join(element) for element in itertools.product(c, repeat=k)]
    ckps = [element for element in itertools.product(mers, repeat=2) if element[0] <= revcomp(element[1])]
    return(ckps)

def save_ckps(ckps, outfile):
    with open(outfile, 'w') as f:
        for ckp in ckps:
            f.write('\t'.join([str(mer2bin(ckp[0])),
                               str(mer2bin(ckp[1])),
                               str(mer2bin(revcomp(ckp[1]))), 
                               str(mer2bin(revcomp(ckp[0]))),
                               ckp[0], 
                               ckp[1], 
                               revcomp(ckp[1]), 
                               revcomp(ckp[0])]) + '\n')
        
def canonical_kmer_pair_main():
    parser = argparse.ArgumentParser(description='canonical k-mer pairs')

    parser.add_argument('--version', action='version', version='%(prog)s 2016-05-25')
    parser.add_argument('-o', metavar = 'o', default = None, required = True,
                        help = 'output file')
    parser.add_argument('-k', metavar = 'k', type=int, required = True,
                        help = 'k')
    parser.add_argument('-e', metavar = 'e', 
                        help = 'eronous k-mer')    
    args = parser.parse_args()

    ckps = canonical_kmer_pairs(args.k, eliminate = args.e)
    save_ckps(ckps, args.o)

if __name__ == '__main__':
    canonical_kmer_pair_main()
