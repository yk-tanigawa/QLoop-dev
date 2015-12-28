import sys, string;

def enumerate_substr(sequence, k):
    substr = set([]);
    for i in range(0, len(sequence) - k + 1):
        substr.add(sequence[i: i + k]);
    return substr;

def replace_N(str_set):
    expanded_set = set([]);
    for x in str_set:
        npos = x.find('N');
        if(npos < 0):
            expanded_set.add(x);
        else:
            rec_set = set([]);
            for n in ['A', 'C', 'G', 'T']:
                rec_set.add(x[:npos] + n + x[npos + 1:]);
            expanded_set.update(replace_N(rec_set));
    return expanded_set;

def prep_seed_set(targets, k):
    substr = set([]);
    for target in targets:
        if(len(target) >= k):
            substr.update(enumerate_substr(target, k));
    return replace_N(substr);

def get_stdin():
    while True:
        try:
            yield ''.join(raw_input());
        except EOFError:
            break;

def seed_match(seedset):
    for line in get_stdin():
        kmerpair = line.split('\t');
        seedmatch = [str(kmer in seedset) for kmer in kmerpair[5:9]]
        print '\t'.join(kmerpair + seedmatch);
        
def revcomp(sequence):
    return string.translate(sequence[::-1], string.maketrans("ACGT", "TGCA"))

def main(argv):
    ctcf = "CCGCGNGGNGGCAG";
    targets = [ctcf, revcomp(ctcf)];
    seed_set = prep_seed_set(targets, 4);
    seed_match(seed_set);

if __name__ == "__main__":
    main(sys.argv);
