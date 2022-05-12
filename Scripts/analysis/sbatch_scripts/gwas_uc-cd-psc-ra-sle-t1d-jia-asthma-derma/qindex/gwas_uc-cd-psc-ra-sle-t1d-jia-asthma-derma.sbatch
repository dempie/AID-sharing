#!/bin/bash

#SBATCH --job-name=qindex_3f
#SBATCH --partition=cpuq
#SBATCH --array=1-335
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=gwas_10-05-2022-%A-%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pietro.demela@fht.org
#SBATCH --mem=64000


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
Rscript qindex_uc_cd_psc_ra_sle_t1d_jia_asthma_derma.script.R  $SLURM_ARRAY_TASK_ID

