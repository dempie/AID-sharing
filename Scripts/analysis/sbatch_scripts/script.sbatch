#!/bin/bash

#SBATCH --job-name=gwas_twof
#SBATCH --partition=cpuq
#SBATCH --array=1-30
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=gwas_twof-%A-%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pietro.demela@fht.org
#SBATCH --mem-per-cpu=
#SBATCH --mem=

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
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID > $SLURM_ARRAY_TASK_ID.txt
Rscript gwas_two_factor_1.R  $SLURM_ARRAY_TASK_ID

