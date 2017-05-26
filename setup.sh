#!/bin/bash

# This will pull all submodules and install LAMMPS
# This script should be run inside the bmr_simulations directory
# After running this script you should be ready to run a simulation

# Update submodules
git submodule update --init --recursive

# Make a .local/bin/ directory
mkdir -p .local/bin/

# Call another script to install LAMMPS
# This script can be run separately to reinstall LAMMPS
bash reinstall_lammps.sh

# Make a symbolic link
cd ../../.local/bin
ln -s lmp_mpi lammps

# Check if the .local/bin/ directory is not on your path
if [[ $PATH != ?(*:)$PWD?(:*) ]]
then
	# Add the .local/bin/ directory to your path
	echo "Adding $PWD to PATH"
	echo -e '\nexport PATH="'$PWD':$PATH"' >> $HOME/.bash_profile
	source $HOME/.bash_profile
fi





