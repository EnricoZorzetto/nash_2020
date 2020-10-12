#!/bin/bash
#SBATCH -o stats.out
#SBATCH -e stats.err
#SBATCH -c 4
#SBATCH --mem=16G
module load GCC/7.4.0
R CMD BATCH --vanilla '--args test=FALSE dataset="G"' main_nash_cluster.R ../Reports/Rout_$SLURM_ARRAY_TASK_ID
# R CMD BATCH --vanilla '--args test=FALSE dataset="G"' main_nash_cluster.R
# R CMD BATCH --vanilla '--args kfold_cv=FALSE dataset="G"' cluster_stats.R
