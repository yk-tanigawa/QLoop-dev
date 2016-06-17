import argparse, itertools

def revcomp(mer, comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}):
    return("".join([comp[mer[i-1:i]] for i in range(len(mer), 0, -1)]))

def mer2bin(mer, table = {'A':'00', 'C':'01', 'G':'10', 'T':'11'}):
    return(int("".join([table[mer[i:i+1]] for i in range(len(mer))]), 2))

def kmers(k, eliminate = None, c = ['A', 'C', 'G', 'T']):
    if(eliminate != None):
        return(["".join(element) for element in itertools.product(c, repeat=k) if eliminate not in "".join(element)])
    else:
        return(["".join(element) for element in itertools.product(c, repeat=k)])
    
def canonical_kmer_pairs(k, eliminate = None, c = ['A', 'C', 'G', 'T']):
    mers = kmers(k, eliminate)
    ckps = [element for element in itertools.product(mers, repeat=2) if element[0] <= revcomp(element[1])]
    return(ckps)

def save_ckps(ckps, outfile):
    if(outfile == None):
        for ckp in ckps:
            print('\t'.join([str(mer2bin(ckp[0])),
                             str(mer2bin(ckp[1])),
                             str(mer2bin(revcomp(ckp[1]))), 
                             str(mer2bin(revcomp(ckp[0]))),
                             ckp[0], 
                             ckp[1], 
                             revcomp(ckp[1]), 
                             revcomp(ckp[0])]))
    else:
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

def save_kmers(mers, outfile):
    if(outfile == None):
        for mer in mers:
            print('\t'.join([str(mer2bin(mer)), mer]))
    else:
        with open(outfile, 'w') as f:
            for mer in mers:
                f.write('\t'.join([str(mer2bin(mer)), mer]) + '\n')
        
def canonical_kmer_pair(k, out, eliminate, pair):
    if(pair):
        ckps = canonical_kmer_pairs(k, eliminate)
        save_ckps(ckps, out)
    else:
        mers = kmers(k, eliminate)
        save_kmers(mers, out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='canonical k-mer pairs')

    parser.add_argument('--version', action='version', version='%(prog)s 2016-06-17')
    parser.add_argument('-o', metavar = 'o', default = None,
                        help = 'output file')
    parser.add_argument('-k', metavar = 'k', type=int, required = True,
                        help = 'k')
    parser.add_argument('-e', metavar = 'e', 
                        help = 'eronous k-mer')    
    parser.add_argument('-p', action='store_true',
                        help = 'k-mer pair (default: false)')
    args = parser.parse_args()

    canonical_kmer_pair(k = args.k,
                        out = args.o,
                        eliminate = args.e,
                        pair = args.p)
