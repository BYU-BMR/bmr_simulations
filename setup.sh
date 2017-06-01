#!/bin/bash

# This will pull all submodules and install LAMMPS
# This script should be run inside the bmr_simulations directory
# After running this script you should be ready to run a simulation

# Update submodules
git submodule update --init --recursive

# Make a .local/bin/ directory
mkdir -p .local/bin/

# Install lammps
cd lammps/src/
make yes-misc
make yes-molecule
make yes-rigid
make yes-user-sph
make yes-user-byu
nice make -j8 mpi

# This deletes the files that were moved into the /src folder
# by the command make yes-user-byu.  Don't worry, this doesn't delete the files 
# out of the USER-BYU directory, running this command just makes it so that
# git doesn't add the files to the repository twice
make no-user-byu

# Copy the executable
cp lmp_mpi ../../.local/bin


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





