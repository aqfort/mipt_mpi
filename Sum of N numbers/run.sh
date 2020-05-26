#!/bin/bash

ARG=101
OUT="task2"

mpic++ -o $OUT $OUT.cpp

while [ ! -e $OUT ]; do
    sleep 1;
done

mpiexec -n 3 ./$OUT ${ARG}

rm $OUT
