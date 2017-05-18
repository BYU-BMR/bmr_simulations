#!/bin/bash


# Some tests showed that using only 6 CPUS all on the same socket was
# faster than using 12 CPUS on the same node.



#SBATCH --mem-per-cpu=2G -t2:00:00


# To submit a bunch of files with similar names:
# for file in `ls in.dpp_pressure_change_0.0089_ea140_ec_100_*` ; do sbatch submit.sh $file ; done

# exit the script if anything has an error
set -e
MACHINEFILE=`/fslapps/fslutils/generate_pbs_nodefile`

#NPUT="$(pwd)/$1"
export OMP_NUM_THREADS=1


constname=2.00E+03_2.72E+03_0_0_1_sh_r_1000_t_0.0005_titlt2_ec_2.00E+02_dc_

dci=1.40
dcf=1.8
difference=0.05

count=`echo "($dcf - $dci)/0.05+1" | bc`
echo " Need  $count  loops"
nn=0;
result=1


for ((a=1; a <= $count ; a++))
do
	dc_need=`echo "($dci+$nn*$difference)" | bc`
	echo $dc_need
	var1=$(awk 'BEGIN{ print "'$dc_need'"<"'$result'" }')  

	if [ "$var1" -eq 1 ];then
	dc_need="0"$dc_need
 	fi

	echo $dc_need
	foldername=$constname$dc_need
 	echo "folder name to go into for data gathering is $foldername"

	cp lammps_data_gather.cpp  $foldername
	cp total_data.txt $foldername

	cd $foldername

	mv  lammps_data_gather.cpp viscosity
	mv total_data.txt  viscosity
	cd   viscosity

	echo "Now in the viscosity folder"
	echo "Gathering data"
		icpc -O3 -std=c++11 lammps_data_gather.cpp -o create_lammps_data_gather

		time ./create_lammps_data_gather	
	echo "End of gathering data"

	let nn=nn+1
	
	mv total_data.txt ..
	cd ..
	mv total_data.txt ..
	cd ..
	echo "Now in the main folder "


done



