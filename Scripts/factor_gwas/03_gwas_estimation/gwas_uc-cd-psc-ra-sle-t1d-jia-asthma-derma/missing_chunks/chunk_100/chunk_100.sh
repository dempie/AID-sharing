#!/bin/bash

#SBATCH --job-name=gw_100
#SBATCH --partition=cpuq
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=gwas_chunk_100_1-10-%A-%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pietro.demela@fht.org
#SBATCH --mem=10000


echo '----------------------------'
echo ' JOB ID: '$SLURM_ARRAY_JOB_ID
echo ' CURRENT TASK ID: '$SLURM_JOB_ID
echo ' CURRENT TASK NUMBER: '$SLURM_ARRAY_TASK_ID
echo '----------------------------'
echo ' MIN TASK ID: '$SLURM_ARRAY_TASK_MIN
echo ' MAX TASK ID: '$SLURM_ARRAY_TASK_MAX
echo ' TOTAL NUMBER OF TASKS: '$SLURM_ARRAY_TASK_COUNT
echo '----------------------------'



module load R/4.1.0
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID 
Rscript chunk_100.R $SLURM_ARRAY_TASK_ID

