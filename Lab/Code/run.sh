#!/bin/bash

OUT="lab"

mpic++ -o $OUT $OUT.cpp

while [ ! -e $OUT ]; do
    sleep 1;
done

mpiexec -np 3 ./$OUT

rm $OUT
