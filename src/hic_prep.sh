#!/bin/sh


CHR="chr21"
RES="1k"
MIN="10k"
MAX="100k"
NORM="KR"
EXP=${NORM}

DATA_DIR="/home/yt/data/GM12878_combined/1kb_resolution_intrachromosomal/${CHR}/MAPQGE30"
BED="/data/yt/Encode/wgEncodeAwgDnaseMasterSites.bed"

python3 ./hic_prep.py \
    -i ${DATA_DIR}/${CHR}_${RES}b.RAWobserved \
    -o ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}.log.norm \
    -c ${CHR} \
    -r ${RES} \
    -m ${MIN} \
    -M ${MAX} \
    -b ${BED} \
    -n ${DATA_DIR}/${CHR}_${RES}b.${NORM}norm \
    -e ${DATA_DIR}/${CHR}_${RES}b.${EXP}expected \
    --log \
    --distnorm
