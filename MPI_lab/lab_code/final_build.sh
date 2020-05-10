#!/bin/bash

mpic++ final.cpp

JOBID=$(qsub final_job.sh)
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

rm -rf a.out
mv $OUT result.txt
rm $ERR

echo "_____ _____ _____ _____ _____"
cat result.txt
echo "_____ _____ _____ _____ _____"
