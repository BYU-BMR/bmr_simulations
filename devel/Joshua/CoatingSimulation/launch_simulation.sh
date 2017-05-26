#!/bin/bash
#
RUN_NAME="$1"
FILENAME=coating
#
echo making directory $RUN_NAME
mkdir $RUN_NAME
echo copying coating.in into $RUN_NAME
cp coating.in $RUN_NAME
echo copying coating.py into $RUN_NAME
cp coating.py $RUN_NAME
echo copying submit.sh into $RUN_NAME 
cp submit.sh $RUN_NAME
#
echo changing directory to $RUN_NAME
cd $RUN_NAME
echo submitting the simulation
sbatch submit.sh