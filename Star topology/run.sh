#!/bin/bash

ARG=(0 1 2 3 4)
OUT="task3"

mpic++ -o $OUT $OUT.cpp

while [ ! -e $OUT ]; do
    sleep 1;
done

mpiexec -n 5 ./$OUT ${ARG[@]}

rm $OUT
