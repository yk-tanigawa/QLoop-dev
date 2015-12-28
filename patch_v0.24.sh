#!/bin/sh

# this is a patch script to generate R plots from 
# results with v0.23

results=$0
results_filtered="${results}.filtered"

cat ${results} | grep -v 'GATC' | \
    python ./ctcf_match.py > ${results_filtered}
    
${DIR}/plots.sh ${results_filtered}
