#!/bin/bash

ARG="0 1 2 3 4 5"

mpic++ task3.cpp

JOBID=$(qsub -v ARG=$ARG job3.sh)

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
