#!/bin/bash

#SBATCH --job-name=gwas_194
#SBATCH --partition=cpuq
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=chunk194_gwasuc_cd_psc_ra_sle_t1d_jia_asthma_derma-%A-%a.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pietro.demela@fht.org
#SBATCH --mem=32000


module load R/4.1.0
Rscript gwas_chunk194_uc_cd_psc_ra_sle_t1d_jia_asthma_derma.script.R  

