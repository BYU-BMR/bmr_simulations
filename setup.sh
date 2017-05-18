#!/bin/bash

# This will pull all submodules and install LAMMPS
# This script should be run inside the bmr_simulations directory
# After running this script you should be ready to run a simulation

# Update submodules
git submodule update --init --recursive

# Make a .local/bin/ directory
mkdir -p .local/bin/

# Build LAMMPS
cd lammps/src/
make yes-misc
make yes-molecule
make yes-rigid
make yes-user-sph
make yes-user-byu
nice make -j8 mpi

# Copy the executable
cp lmp_mpi ../../.local/bin

# Make a symbolic link
cd ../../.local/bin
ln -s lmp_mpi lammpz

# Check if the .local/bin/ directory is not on your path
if [[ $PATH != ?(*:)$PWD/.local/bin?(:*) ]]
then
	# Add the .local/bin/ directory to your path
	echo "Adding $PWD/.local/bin to PATH"
	echo -e '\nexport PATH="'$PWD'/.local/bin:$PATH"' >> $HOME/.bash_profile
	source $HOME/.bash_profile
else





