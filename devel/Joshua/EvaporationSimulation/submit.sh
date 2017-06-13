#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=12   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=4096M   # memory per CPU core
#SBATCH -J "quick_coat"   # job name
#MACHINEFILE=`/fslapps/fslutils/generate_pbs_nodefile`

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=1

module load python/3/4

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
time python coating.py && mpirun -np $SLURM_NTASKS lammps -in solvent.in