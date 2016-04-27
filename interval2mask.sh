#!/bin/sh

BEDGZ="/data/yt/Encode/wgEncodeAwgDnaseMasterSites.bed.gz"
CHR=21
RES=1000
MARGIN=1000

gzip -dc ${BEDGZ} \
    | grep chr${CHR} \
    | cut -f 2,3 \
    | ./interval2mask.pl --res=${RES} --margin=${MARGIN}
