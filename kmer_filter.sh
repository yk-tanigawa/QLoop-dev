#!/bin/sh

# Mail at end/on suspension
#$ -m es

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q 24core.q
#$ -pe smp 24

set -eu

CHR="chr21"
RES="1k"
MIN="10k"
MAX="20k"
NORM="SQRTVC"
MARGIN=0
EXP=${NORM}
k=7
ELIMINATE="GATC"
iter1=2000
iter2=100
acc=0.1

PERCENTILE_CIRCOS=0.95

VERSION="v0.60"
PROG_DIR="/work2/yt/QLoop-dev"
DATA_DIR="/data/yt/GM12878_combined/1kb_resolution_intrachromosomal/${CHR}/MAPQGE30"
RESULTS_DIR="/data/yt/QLoop-dev/${VERSION}"
FASTA="/data/yt/hg19/chromosomes/${CHR}.fa"
KMER="/data/yt/QLoop-dev/canonical_kmer_pair"
BED="/data/yt/Encode/wgEncodeAwgDnaseMasterSites.bed"

##PRI="/data/yt/QLoop-dev/v0.57/v057_chr21_SQRTVC_mar0_m10M20_acc01.sh.e424194.res"

# output results file
RES_FILE="${RESULTS_DIR}/${VERSION}_${CHR}_${NORM}_mar${MARGIN}_m${MIN}M${MAX}_nu${acc}_iter${iter1}.res"
CMP="${RES_FILE}.self"


##################################################
# make a copy of this script
##################################################
#cat $0 > "${SGE_O_WORKDIR}/${JOB_NAME}.s${JOB_ID}"

##################################################
# prreparation step
##################################################
if [ ! -d ${RESULTS_DIR} ]; then
    # create a directory for the results file
    mkdir -p ${RESULTS_DIR} 
fi

if [ ! -e ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm ]; then
    # pre-processing of a Hi-C matrix
    python3 ${PROG_DIR}/src/hic_prep.py \
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
fi

##################################################
# git:
#  switch version and rebuild
##################################################
cd ${PROG_DIR}

#git checkout ${VERSION}
#git log --oneline --graph --decorate -n3

make clean
make kmer_filter

${PROG_DIR}/kmer_filter -v

##################################################
# perform twin Boost
#  (current implementation is L2 Boosting)
##################################################
${PROG_DIR}/kmer_filter \
    -k ${k} \
    --res 1000 \
    --margin ${MARGIN} \
    --iter1 ${iter1} \
    --iter2 ${iter2} \
    --acc ${acc} \
    --fasta ${FASTA} \
    --hic ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm \
    --kmer ${KMER}/k${k}.e${ELIMINATE}.kmers \
    --out ${RES_FILE} # --pri ${PRI} 

git checkout dev
