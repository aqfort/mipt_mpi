#!/bin/bash

OUT="task1"

mpic++ -o $OUT $OUT.cpp

while [ ! -e $OUT ]; do
    sleep 1;
done

mpiexec -n 10 ./$OUT

rm $OUT
