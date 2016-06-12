#!/bin/sh

set -eu

CIRCOS_TEMP_DIR=`basename ${CSVFILE}`.tmp
CIRCOS_BIN_DIR="/work2/yt/usr/circos/current/bin"

cat ${CSVFILE} \
    | sed -e 's/^ /- /' \
    | perl -I/home/yt/perl5/lib/perl5/ ${CIRCOS_BIN_DIR}/circos-tableviewer-parse-table -conf /work2/yt/QLoop-dev/src/parse-table.conf \
    | perl -I/home/yt/perl5/lib/perl5/ ${CIRCOS_BIN_DIR}/circos-tableviewer-make-conf -dir ${CIRCOS_TEMP_DIR}

echo "circos table viwer complete!"

cat /work2/yt/QLoop-dev/src/circos.conf \
    | sed -e "s/__DATA__/${CIRCOS_TEMP_DIR}/" \
    | sed -e "s!__IMGDIR__!${OUT_DIR}!" > ${CIRCOS_TEMP_DIR}/circos.conf 

perl -I/home/yt/perl5/lib/perl5/ ${CIRCOS_BIN_DIR}/circos \
    -conf ${CIRCOS_TEMP_DIR}/circos.conf \
    -outputfile $(basename ${CSVFILE}).png \
    | grep created
    
rm -rf ${CIRCOS_TEMP_DIR}
