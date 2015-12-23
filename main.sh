#!/bin/sh

make clean
make main

prog_name="ChromLoopC"
version="v0.20"

DIR="."
DATA_DIR_ROOT="/data/yt"


# parameters
chr=21
k=2
res=1000
min_size=10000
max_size=100000
iteration_num=100
percentile=0.5
norm="KR"
exp=${norm}
genome="GRCh37"

##########################
# input file names
fasta_file="${DATA_DIR_ROOT}/${genome}.ch${chr}.fasta"
hicRaw_dir="${DATA_DIR_ROOT}/GM12878_combined"
kmerFreq_file="./tmp/GRCh37.ch21.r1000.k2.freq"
hic_file="./tmp/chr21.m10000.M100000.VC.VC.hic"
boostOracle_file=".stamps"

# output
output_dir="./${version}/"


${DIR}/main \
    --chr ${chr} \
    --k ${k} \
    --res ${res} \
    --min_size ${min_size} \
    --max_size ${max_size} \
    --iteration_num ${iteration_num} \
    --percentile ${percentile} \
    --norm ${norm} \
    --expected ${exp} \
    --fasta ${fasta_file} \
    --hicRaw ${hicRaw_dir} \
    --kmerFreq ${kmerFreq_file} \
    --hic ${hic_file} \
    --boostOracle ${boostOracle_file} \
    --out ${output_dir} \
    --skipPrep
  