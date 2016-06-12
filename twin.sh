#!/bin/sh

# Mail at end/on suspension
#$ -m es

#$ -S /bin/sh
#$ -cwd
#$ -V
# -q 24core.q
#$ -q blade.q
#$ -pe smp 12

set -eu

CHR="chr21"
RES="1k"
MIN="10k"
MAX="20k"
NORM="SQRTVC"
MARGIN=0
EXP=${NORM}
k=5
ELIMINATE="GATC"
iter1=401
iter2=100
acc=0.1

PERCENTILE_CIRCOS=0.95

VERSION="v0.58"
PROG_DIR="/work2/yt/QLoop-dev"
DATA_DIR="/data/yt/GM12878_combined/1kb_resolution_intrachromosomal/${CHR}/MAPQGE30"
RESULTS_DIR="/data/yt/QLoop-dev/${VERSION}"
FASTA="/data/yt/hg19/chromosomes/${CHR}.fa"
KMER="/data/yt/QLoop-dev/canonical_kmer_pair"
BED="/data/yt/Encode/wgEncodeAwgDnaseMasterSites.bed"

PRI="/data/yt/QLoop-dev/v0.57/v057_chr21_SQRTVC_mar0_m10M20_acc01.sh.e424194.res"

# output results file
RES_FILE="${RESULTS_DIR}/${VERSION}_${CHR}_${NORM}_mar${MARGIN}_m${MIN}M${MAX}_nu${acc}_iter${iter1}.res"
CMP="${RES_FILE}.self"


##################################################
# make a copy of this script
##################################################
cat $0 > "${SGE_O_WORKDIR}/${JOB_NAME}.s${JOB_ID}"

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
git checkout dev

git log --oneline --graph --decorate -n3

make clean
make twin

${PROG_DIR}/twin -v

##################################################
# perform twin Boost
#  (current implementation is L2 Boosting)
##################################################
${PROG_DIR}/twin \
    -k ${k} \
    --res 1000 \
    --margin ${MARGIN} \
    --iter1 ${iter1} \
    --iter2 ${iter2} \
    --acc ${acc} \
    --fasta ${FASTA} \
    --hic ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm \
    --kmer ${KMER}/k${k}.e${ELIMINATE}.ckp \
    --pri ${PRI} \
    --out ${RES_FILE}

##################################################
# plot summary of the results
##################################################
python3 ${PROG_DIR}/src/plt_res.py \
    -i ${RES_FILE} \
    -o ${RES_FILE}.png

##################################################
# circos plot
##################################################
python3 ${PROG_DIR}/src/res_dump.py \
    -r ${RES_FILE} \
    -c ${KMER}/k${k}.e${ELIMINATE}.ckp \
    -p ${PERCENTILE_CIRCOS} \
    -o ${RES_FILE}.p${PERCENTILE_CIRCOS}.csv

OUT_DIR=${RESULTS_DIR} CSVFILE="${RES_FILE}.p${PERCENTILE_CIRCOS}.csv" ${PROG_DIR}/src/make_circos.sh

##################################################
# git:
#  switch version and rebuild
##################################################
cd ${PROG_DIR}

#git checkout ${VERSION}
git checkout dev

git log --oneline --graph --decorate -n3

make clean
make pred

##################################################
# self prediction
##################################################

${PROG_DIR}/pred \
    -k ${k} \
    --res 1000 \
    --margin ${MARGIN} \
    --fasta ${FASTA} \
    --hic ${DATA_DIR}/${CHR}_${RES}b_m${MIN}_M${MAX}_${NORM}_${EXP}_mar${MARGIN}.log.norm \
    --kmer ${KMER}/k${k}.e${ELIMINATE}.ckp \
    --pri ${RES_FILE} \
    --out ${CMP}

python3 ${PROG_DIR}/src/scatter_plot.py \
    -i ${CMP}.cmp

git checkout dev
