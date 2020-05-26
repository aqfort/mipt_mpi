#!/bin/bash

ARG=101

mpic++ task2.cpp

JOBID=$(qsub -v ARG=$ARG job2.sh)

echo "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░"

echo Job ID: $JOBID

TEMP=$(echo $JOBID | cut -f1 -d.)

OUT="000.o$TEMP"
ERR="000.e$TEMP"

while [ ! -e $OUT ]; do
    sleep 1;
done

while [ ! -e $ERR ]; do
    sleep 1;
done

mv $OUT result.txt
rm $ERR

echo "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░"
cat result.txt
echo "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░"
