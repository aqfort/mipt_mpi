#!/bin/bash

ARG=10000000
OUT="task4"

mpic++ -o $OUT $OUT.cpp

while [ ! -e $OUT ]; do
    sleep 1;
done

mpiexec -n 3 ./$OUT ${ARG}

rm $OUT
