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
FASTA="/data/yt/hg19/chromosomes/${CHR}.fa.gz"
BED="/data/yt/Encode/wgEncodeAwgDnaseMasterSites.bed"

python3 /work2/yt/QLoop-dev/src/hic_prep.py \
    -i ${DATA_DIR}/${CHR}_${RES}b.RAWobserved \
    -o ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm \
    -c ${CHR} \
    -r ${RES} \
    --margin ${MARGIN} \
    -m ${MIN} \
    -M ${MAX} \
    -f ${FASTA} \
    -b ${BED} \
    -n ${DATA_DIR}/${CHR}_${RES}b.${NORM}norm \
    -e ${DATA_DIR}/${CHR}_${RES}b.${EXP}expected \
    --log \
    --distnorm 
