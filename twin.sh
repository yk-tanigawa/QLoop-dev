#!/bin/sh

# Mail at end/on suspension
#$ -m es

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q

CHR="chr21"
RES="1k"
MIN="10k"
MAX="100k"
NORM="KR"
MARGIN=500
EXP=${NORM}

DATA_DIR="/data/yt/GM12878_combined/1kb_resolution_intrachromosomal/${CHR}/MAPQGE30"
FASTA="/data/yt/hg19/chromosomes/${CHR}.fa"

./twin \
    -k 5 \
    --res 1000 \
    --margin ${MARGIN} \
    --iter1 10000 \
    --iter2 1000 \
    --eliminate GATC \
    --fasta ${FASTA} \
    --hic ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm \
    --out ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm 
