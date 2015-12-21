#!/bin/sh
cat $0

DIR="/work2/yt/ChromLoopC"
cd ${DIR}

prog_version="20151221"

percentile=0.5
norm="KR"
T=100
k=4
m=10000
M=100000
exp=${norm}
res=1000
chr=21

BoostMode=""

HiCdataFile=${DIR}/data/chr${chr}.m${m}.M${M}.${norm}.${exp}.dat
BoostFile=${DIR}/data/${prog_version}/chr${chr}.m${m}.M${M}.${norm}.${exp}.dat.k${k}.t${threshold}.T${T}${BoostMode}.stamps
kmerPairFile=${BoostFile}.kmerpair
QPfileP=${kmerPairFile}.P${T}.dat
QPfileq=${kmerPairFile}.q${T}.dat
QPfileOut=${kmerPairFile}.qpout

git checkout v${prog_version}
make clean
make QPprep kmerPairBoost

threshold=`./hist.sh ${HiCdataFile}|awk -v percentile="${percentile}" '{if($1 == percentile) print $2}'`

${DIR}/kmerPairBoost \
    -f ${DIR}/data/k${k}.res${res}.chr${chr}.dat \
    -H ${HiCdataFile} \
    -k${k} \
    -r${res} \
    -T${T} \
    -t${threshold} \
    -o ${DIR}/data/${prog_version}

cat ${BoostFile} |awk '{print $4}' > ${kmerPairFile}

${DIR}/QPprep \
    --frequency ${DIR}/data/k${k}.res${res}.chr${chr}.dat \
    --hic ${HiCdataFile} \
    --kmerpair ${KmerPairFile} \
    -k ${k} \
    --res ${res} \
    --out ${DIR}/data/${prog_version}

/work2/yt/QuadProg-example/QPwithFile ${QPfileP} ${QPfileq} ${T} > ${QPfileOut}

#################################
