#!/bin/bash
#
#SBATCH -p serial_requeue # partition (queue)
#SBATCH -c 6 # number of cores
#SBATCH --mem 96GB 
#SBATCH -t 4-08:00 # time (D-HH:MM)
#SBATCH --array=2000-2016

singularity exec --cleanenv --env R_LIBS_USER=../RELEASE_3_14 ../RELEASE_3_14.sif Rscript code/sim_original.R -y $SLURM_ARRAY_TASK_ID -m vanilla2 -w 0.75
 
