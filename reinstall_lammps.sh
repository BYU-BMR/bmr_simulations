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
