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
NORM="SQRTVC"
MARGIN=0
EXP=${NORM}
k=5
ELIMINATE="GATC"
f_norm="L1"

version="v0.61"
PROG_DIR="/work2/yt/QLoop-dev"
DATA_DIR="/data/yt/GM12878_combined/1kb_resolution_intrachromosomal/${CHR}/MAPQGE30"
FASTA="/data/yt/hg19/chromosomes/${CHR}.fa"
KMER="/data/yt/QLoop-dev/canonical_kmer_pair"
PRI="/data/yt/QLoop-dev/v0.54/v054N400.e413446.res"

cd ${PROG_DIR}

#git checkout ${version}

#git log --oneline --graph --decorate -n3

#make clean
#make twin

${PROG_DIR}/pred -v

${PROG_DIR}/pred \
    -k ${k} \
    --res 1000 \
    --margin ${MARGIN} \
    --fasta ${FASTA} \
    --hic ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm \
    --kmer ${KMER}/k${k}.e${ELIMINATE}.ckp \
    --pri ${PRI} \
    --f_norm ${f_norm} \
    --out ${PRI}.${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm

#git checkout master
