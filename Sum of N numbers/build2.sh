#!/bin/bash
mpic++ 2prog.cpp
qsub -v ARG=101 job2.sh
