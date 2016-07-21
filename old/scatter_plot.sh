#!/bin/sh

# Mail at end/on suspension
#$ -m es

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q

PROG_DIR="/work2/yt/QLoop-dev"
CMP_FILE="/data/yt/QLoop-dev/v0.54/v054N400.e413446.res.chr21_1kb_m10k_M100k_KR_KR_mar500.log.norm.cmp"
cat $0

python3 ${PROG_DIR}/src/scatter_plot.py -i ${CMP_FILE} -o ${CMP_FILE}.png

#################################
