#!/bin/sh

# Mail at end/on suspension
#$ -m es

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q

cat $0

k=5
e="GATC"

python3 /work2/yt/QLoop-dev/src/canonical_kmer_pair.py \
    -o /data/yt/QLoop-dev/canonical_kmer_pair/k${k}.e${e}.ckp \
    -k ${k} \
    -e ${e}

#################################
