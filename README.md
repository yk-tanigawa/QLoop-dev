# QLoop-dev
QLoop: a novel quantitative method to analyze Hi-C matrix to infer critical factors in chromatin looping

# Usage
$./twin \
       -k k \
       --res r \
       [--margin M] \
       --iter1 n \
       --iter2 m \
       [--acc a] \
       --fasta f \
       --hic H \
       --kmer c \
       --out o \
       [--pri p] \
       [--sec s] \
       [--verbose V] \
       --thread_num t

- k : kmer-length
- r : resolution
- M : margin to count k-mer frequency
- n : iteration num. in the first round of twin boosting
- m : iteration num. in the second round of twin boosting
- a : acceleration paremeter in L2 Boosting (0 < a <= 1.0)
- f : fasta file (now only supports unzipped file as of v0.51)
- H : pre-processed Hi-C file
      you can pre-process Hi-C raw file with src/hic_prep.py
- c : canonical k-mer pair file
- o : output file name (unsupported as of v0.51)
- p : saved results of the first round of twin boosting (unsupported as of v0.51)
- s : saved results of the second round of twin boosting (unsupported as of v0.51)
- V : verbose level (unsupported as of v0.51)
- t : thread num
