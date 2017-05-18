#!/bin/bash


# Some tests showed that using only 6 CPUS all on the same socket was
# faster than using 12 CPUS on the same node.



#SBATCH --mem-per-cpu=2G -t48:00:00 


# To submit a bunch of files with similar names:
# for file in `ls in.dpp_pressure_change_0.0089_ea140_ec_100_*` ; do sbatch submit.sh $file ; done

# exit the script if anything has an error
set -e
MACHINEFILE=`/fslapps/fslutils/generate_pbs_nodefile`

#NPUT="$(pwd)/$1"
module load lammps
module load matlab
export OMP_NUM_THREADS=1 

for i in {2..2}
do
   	echo "series $i times"

	foldername="$i"

	if [ ! -d "$foldername" ]; then

	mkdir $foldername

	fi

	mv number.txt $foldername
#	cp number.txt $foldername
        cp lj_potential.txt $foldername
        cp diameter_Active.txt $foldername
	cp sbatch_and_potential_creator_active_solid.cpp $foldername
	cd $foldername

	
		icpc -O3 -std=c++11 sbatch_and_potential_creator_active_solid.cpp -o create_sbatch_and_potential

		time ./create_sbatch_and_potential	


	echo "finish series $i"
	
	mv number.txt ..

	#mv slurry.sh ..
	#mv potential.mod ..
	#mv potential1.mod ..
	#mv lammps_input_parameters.txt ..
	

        sbatch slurry.sh


	cd ..

	sleep 1m



done




