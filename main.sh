#!/bin/sh

make clean
make main

DIR="."

chr=21
k=2
res=1000
min_size=10000
max_size=100000
iteration_num=100
percentile=0.5
norm="KR"
exp=${norm}


${DIR}/main \
    --chr ${chr} \
    --k ${k} \
    --res ${res} \
    --min_size ${min_size} \
    --max_size ${max_size} \
    --iteration_num ${iteration_num} \
    --percentile ${percentile} \
    --norm ${norm} \
    --expected ${exp}
  