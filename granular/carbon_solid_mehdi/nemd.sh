#!/bin/bash


# Some tests showed that using only 6 CPUS all on the same socket was
# faster than using 12 CPUS on the same node.
#SBATCH --nodes=1 --ntasks=24
#SBATCH --mem-per-cpu=2G -t96:00:00 


# To submit a bunch of files with similar names:
# for file in `ls in.dpp_pressure_change_0.0089_ea140_ec_100_*` ; do sbatch submit.sh $file ; done

# exit the script if anything has an error
set -e


#NPUT="$(pwd)/$1"

INPUT=in.nemd 
foldername=viscosity

SUCC="success.txt"
MACHINEFILE=`/fslapps/fslutils/generate_pbs_nodefile`


module load lammps
module load matlab


export OMP_NUM_THREADS=1

if [ ! -d "$foldername" ]; then

mkdir $foldername

fi


cp $INPUT $foldername

cp restart.viscosity_elasticity_calculation $foldername

cp potential.mod  $foldername
cp potential1.mod  $foldername


cd $foldername

	
###while loop


#while [ ! -f "$SUCC" ]
#do
	

	echo "initialing lammps"

	mpirun -np $SLURM_NTASKS -machinefile $MACHINEFILE $(which lammps) <$INPUT

	
	echo "finish lammps"
	

#done


cd ..
scp slurm-$SLURM_JOB_ID.out $foldername

