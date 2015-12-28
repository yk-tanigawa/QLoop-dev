#!/bin/sh

R --vanilla --slave \
    --args $1 \
    < histo.R 
