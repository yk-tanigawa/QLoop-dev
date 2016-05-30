#!/bin/sh

# Mail at end/on suspension
#$ -m es

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q 24core.q
#$ -pe smp 24

CHR="chr21"
RES="1k"
MIN="10k"
MAX="100k"
NORM="KR"
MARGIN=500
EXP=${NORM}
k=5
ELIMINATE="GATC"
iter1=1000
iter2=100
acc=1.0

version="v0.52"
PROG_DIR="/work2/yt/QLoop-dev"
DATA_DIR="/data/yt/GM12878_combined/1kb_resolution_intrachromosomal/${CHR}/MAPQGE30"
FASTA="/data/yt/hg19/chromosomes/${CHR}.fa"
KMER="/data/yt/QLoop-dev/canonical_kmer_pair"


cd ${PROG_DIR}

git checkout ${version}

git log --oneline --graph --decorate -n3

make clean
make twin

${PROG_DIR}/twin -v

echo ${PROG_DIR}/twin \
    -k ${k} \
    --res 1000 \
    --margin ${MARGIN} \
    --iter1 ${iter1} \
    --iter2 ${iter2} \
    --acc ${acc} \
    --fasta ${FASTA} \
    --hic ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm \
    --kmer ${KMER}/k${k}.e${ELIMINATE}.ckp \
    --out ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm.k${k}

git checkout master