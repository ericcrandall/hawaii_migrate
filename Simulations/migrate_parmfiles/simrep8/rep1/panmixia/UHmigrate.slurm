#!/bin/bash
#SBATCH -J migrateUH # Name for your job
#SBATCH -n 1 # Number of cores requested, a node has 20 cores in the all.q partition
#SBATCH -N 1 # Number of nodes, Current Maximum is 10 nodes in the all.q partition - requests over this maximum will not run
#SBATCH -t 43200 # Runtime in minutes. The Maximum runtime currently is 72 hours, 4320 minutes - requests over that time will not run
#SBATCH -p longtest.q # Partition to submit to the standard compute node partition(all.q) or the large memory node partition(lm.q)
#SBATCH -o migrateUH.out # Standard out goes to this file
#SBATCH -e migrateUH.err # Standard err goes to this file
#SBATCH --mail-user eric.d.crandall@gmail.com # this is the email you wish to be notified at
#SBATCH --mail-type ALL # this specifies what events you should get an email about ALL will alert you of job beginning,completion, failure etc
#source ~/.bash_profile #if you want to use modules or need environment variables use this if your shell is bash to load those
../../../../migrate-3.6.6/migrate-n parmfile #put the commands you want to run