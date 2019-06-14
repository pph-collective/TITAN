#!/bin/bash

for i in 1 3 5 7 13
do
    ./subTitan.sh calTrt.py -N 55000 -j memory_$i -m $i
done
