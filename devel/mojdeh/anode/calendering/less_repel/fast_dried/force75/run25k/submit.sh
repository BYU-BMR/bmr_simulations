#!/bin/bash

#SBATCH --time=72:00:00   # walltime
#SBATCH --ntasks=12   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem=4096M   # memory per CPU core
#SBATCH -J "75-25kless_repel"   # job name
#MACHINEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
#SBATCH --mail-user=mojdeh_n87@yahoo.com   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Compatibility variables for PBS. Delete if not needed.
export PBS_NODEFILE=`/fslapps/fslutils/generate_pbs_nodefile`
export PBS_JOBID=$SLURM_JOB_ID
export PBS_O_WORKDIR="$SLURM_SUBMIT_DIR"
export PBS_QUEUE=batch

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=1

module load python/3/4

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# Load appropiate modules
# module purge
# module load lammps/21Oct14
# module load lammps/2015-05-15_openmpi-1.8.5_gnu-5.1.0
# module load compiler_intel/15.0.0
# module load mkl/11.2.0
# module load mpi/openmpi-1.8.1_gnu-4.9.0
# module load compiler_gnu/4.9.0
# module load gdb/7.9.1
# module load matlab
# module load python/2/7
# module load compiler_gnu/4.9
# module load mpi/openmpi-1.8_gnu-4.9


# Compile sbatch_and_potential_creator.cpp and run executable
#python generateWallAtoms_new.py
# python tabular_generator.V1.py
# python small_box_new_viscosity-expansion.py

# Run in.granular
#mpirun -np 16 lammps -in in.granular
#lammps -in in.granular
# mpirun -np $SLURM_NTASKS /fslgroup/fslg_bmr_wheeler/lammps-30Jul16/src/lmp_mpi -in in.drop
time python calendering.py && mpirun -np $SLURM_NTASKS lammps -in calendering.in