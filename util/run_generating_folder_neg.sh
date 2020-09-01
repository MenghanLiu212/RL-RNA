#!/bin/bash

####SBATCH -p serial                   # Send this job to the serial partition
#SBATCH -n 4                        # number of cores
####SBATCH -q wildfire                # send to wildfire queue
#SBATCH -t 5-12:50                 # wall time (D-HH:MM)
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=mliu126@asu.edu # send-to address

module load viennarna/2.4.14
module load python/2.7.15

cd /home/mliu126/RNA/RL-RNA/Gnerating_neg_label_folder/ENTRNA-master_new_MH/ENTRNA-master/util
python generating_folder_for_training.py
