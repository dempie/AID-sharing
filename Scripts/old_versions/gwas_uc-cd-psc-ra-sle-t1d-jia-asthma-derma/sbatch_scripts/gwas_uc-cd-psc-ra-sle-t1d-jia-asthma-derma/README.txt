sbatch_scripts/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma

This folder contains the script that have been used to send the job and job arrays to run the userGWAS with genomiSEM.

a) chunks/ contains a copy of the sbatch and R.script to run the userGWAS with the model specified. 
b) qindex/ contains a copy of the sbatch and R.script to run the userGWAS with the model for the Q index calculation

c) chunk_12/ contains a copy of the sbatch and R.script to run the userGWAS only on chunk_12 
c) chunk_25/ contains a copy of the sbatch and R.script to run the userGWAS only on chunk_25
c) chunk_87/ contains a copy of the sbatch and R.script to run the userGWAS only on chunk_87 
c) chunk_131/ contains a copy of the sbatch and R.script to run the userGWAS only on chunk_131
c) chunk_194/ contains a copy of the sbatch and R.script to run the userGWAS only on chunk_194
c) chunk_203/ contains a copy of the sbatch and R.script to run the userGWAS only on chunk_203 

The folders for the chunks 12,25,87,131,194,203 were necessary because this chunks produced an error during the esitmation in the job-array sent with the scripts in chunks/. The errors were different, from the log file, the erro for each of the chunk was: 

#chunk 12 : task 144 failed - "system is computationally singular: reciprocal condition number = 3.33818e-38"
#chunk25 : Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 87 : Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 131: task 211 failed - "system is computationally singular: reciprocal condition number = 1.50212e-24"
#chunk 194: Error in unserialize(socklist[[n]]) : error reading from connection. for this just sent the job again
#chunk 203:  task 148 failed - "system is computationally singular: reciprocal condition number = 2.86735e-37"


Log files can be found ain /AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks

The reason for the computationally singular error we think is due for issues estimating the model that produces matrixes that cannot be inverted. Thus I decided to split the chunks in subchunks in order to identify a group of SNPs that causes the system to be singular. Results of this splitting of subchunks at /AID_sharing/outputs/gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/05_gwas_ouput/chunks/

I just ignore the subchunks that give an error after the subsplitting of the chunks 12,25,87,131,194,203.

The reason for the error 'Error in unserialize(socklist[[n]]) : error reading from connection' for this just sent the job again is because the hpc, so it's sufficient to re-run again the chunks
