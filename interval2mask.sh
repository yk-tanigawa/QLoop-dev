#!/bin/sh

#$ -m e

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q

BED="/data/yt/Encode/wgEncodeAwgDnaseMasterSites"
BEDGZ="${BED}.bed.gz"
CHR=21
RES=1000
MARGIN=0
OUT="${BED}.chr${CHR}.res${RES}.mar${MARGIN}.interval"

gzip -dc ${BEDGZ} \
    | grep chr${CHR} \
    | cut -f 2,3 \
    | ./interval2mask.pl --res=${RES} --margin=${MARGIN} > ${OUT}
