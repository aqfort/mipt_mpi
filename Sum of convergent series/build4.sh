#!/bin/bash

ARG=10000000

mpic++ 4prog.cpp

# ##### ##### ##### ##### #####
# COMMANDS="job4.sh"
# JOBIDS=""
# for CMD in $COMMANDS; do
#     JOBIDS="$JOBIDS:`qsub $CMD`"
# done
# qsub -W depend=afterok:$JOBIDS postprocessing.sh
# ##### ##### ##### ##### #####

JOBID=$(qsub -v ARG=$ARG job4.sh)
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
