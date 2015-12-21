#!/bin/sh

R --vanilla --slave \
    --args $1 \
    < kmerPairPlot.R 
