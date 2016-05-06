#!/bin/sh

#$ -m bes

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q 24core.q

#$ -l fc=18
#$ -l mf=2G

########################################################################
# This is a wrapper script to run the main program of QLoop-dev
########################################################################

prog_name="QLoop-dev"
version="v0.41"

# parameters
chr=21
k=5
res=1000
min_size=10000
max_size=100000
iteration_num=1000
percentile=0.8
percent=80
norm="KR"
exp=${norm}
DATA_DIR_ROOT="/data/yt"
genome="GRCh37"

DIR="/work2/yt/${prog_name}"

##########################
# input file names
fasta_file="${DATA_DIR_ROOT}/${genome}.ch${chr}.fasta"
hicRaw_dir="${DATA_DIR_ROOT}/GM12878_combined"
kmerFreq_file="${DIR}/tmp/${genome}.ch${chr}.r${res}.k${k}.freq"
hic_file="${DIR}/tmp/chr21.m${min_size}.M${max_size}.${norm}.${exp}.hic"
boostOracle_file="hoge.oracle"
interval_file="/home/yt/data/Encode/wgEncodeAwgDnaseMasterSites.chr21.res1000.mar0.interval"

# output
output_dir="${DATA_DIR_ROOT}/${prog_name}/${version}/"

max_kb=`expr ${max_size} / 1000`
min_kb=`expr ${min_size} / 1000`
res_kb=`expr ${res} / 1000`

histo="${output_dir}chr21.m${min_kb}k.M${max_kb}k.${norm}.${exp}.histo"
stamps="${output_dir}chr21.m${min_kb}k.M${max_kb}k.${norm}.${exp}.k${k}.res${res_kb}k.p${percent}.T${iteration_num}.stamps"
QP_P="${output_dir}chr21.m${min_kb}k.M${max_kb}k.${norm}.${exp}.k${k}.res${res_kb}k.p${percent}.T${iteration_num}.P"
QP_q="${output_dir}chr21.m${min_kb}k.M${max_kb}k.${norm}.${exp}.k${k}.res${res_kb}k.p${percent}.T${iteration_num}.q"
QP_out="${output_dir}chr21.m${min_kb}k.M${max_kb}k.${norm}.${exp}.k${k}.res${res_kb}k.p${percent}.T${iteration_num}.QP"
results="${output_dir}chr21.m${min_kb}k.M${max_kb}k.${norm}.${exp}.k${k}.res${res_kb}k.p${percent}.T${iteration_num}.results"
results_filtered="${output_dir}chr21.m${min_kb}k.M${max_kb}k.${norm}.${exp}.k${k}.res${res_kb}k.p${percent}.T${iteration_num}.results.filtered"

##########################

cat $0
cd ${DIR}

git log --oneline --graph --decorate -n3
git checkout ${version}
git log --oneline --graph --decorate -n3
make clean
make main

if [ ! -e ${output_dir} ]; then mkdir -p ${output_dir}; fi

${DIR}/main -v

${DIR}/main \
    --chr ${chr} \
    --k ${k} \
    --res ${res} \
    --min_size ${min_size} \
    --max_size ${max_size} \
    --iteration_num ${iteration_num} \
    --percentile ${percentile} \
    --norm ${norm} \
    --expected ${exp} \
    --fasta ${fasta_file} \
    --hicRaw ${hicRaw_dir} \
    --kmerFreq ${kmerFreq_file} \
    --hic ${hic_file} \
    --boostOracle ${boostOracle_file} \
    --interval ${interval_file} \
    --out ${output_dir} \
    --skipPrep 

${DIR}/histo.sh ${histo}

/work2/yt/QuadProg-example/QPwithFile \
    ${QP_P} ${QP_q} ${iteration_num} > \
    ${QP_out}

paste ${QP_out} ${stamps} > ${results}

cat ${results} | grep -v 'GATC' | \
    python ./ctcf_match.py ${k} > ${results_filtered}
    
${DIR}/plots.sh ${results_filtered}

git checkout master
