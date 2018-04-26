#! /bin/bash
count=0
for model in */; do

	echo "$model"
	cd "$model"
	#rm migrate.* 
	sbatch UHmigrate.slurm
	count=`expr $count + 1`
	cd ..
	done
	
echo "done submitting $count runs"

