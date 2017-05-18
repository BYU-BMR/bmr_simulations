#!/bin/bash

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

# Add the .local/bin/ directory to your path
echo -e '\nexport PATH="'$PWD'/.local/bin:$PATH"' >> $HOME/.bash_profile
source $HOME/.bash_profile


