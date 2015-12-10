#!/bin/sh

awk '{print $3}' $1 | ./hist_sub
