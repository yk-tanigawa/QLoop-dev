#!/bin/sh

# this is a patch script to generate R plots from 
# results with v0.23

DIR=`pwd`

results="${DIR}/$1"
results_filtered="${results}.filtered"
k="${2}"


cat ${results} | grep -v "GATC" | \
    python ./ctcf_match.py ${k} > ${results_filtered}
    
${DIR}/plots.sh ${results_filtered}
