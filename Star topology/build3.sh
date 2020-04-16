#!/bin/bash
mpic++ 3prog.cpp
qsub -v ARG="0 1 2 3 4 5" job3.sh
