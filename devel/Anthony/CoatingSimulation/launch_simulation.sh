#!/bin/bash
#
RUN_NAME="$1"
FILENAME=test
#
echo making directory $RUN_NAME
mkdir $RUN_NAME
echo copying test.in into $RUN_NAME
cp test.in $RUN_NAME
echo copying test.py into $RUN_NAME
cp test.py $RUN_NAME
echo copying submit.sh into $RUN_NAME 
cp submit.sh $RUN_NAME
#
echo changing directory to $RUN_NAME
cd $RUN_NAME
echo submitting the simulation
sbatch submit.sh