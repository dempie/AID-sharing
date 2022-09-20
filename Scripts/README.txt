Subfolders of Scripts/

a) test/ contains script with tests and script where I was learning how to use the package 
b) summary_stats_info/ contains a script which outputs a table with info on the GWAS used here
c) functions/ contains scripts to store my functions
  
  
d) gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/ in this folder the code for running and assemblying this GWAS
          it was a 3 factor GWAS with the following diseases: 
              ulcerative colitis, (de Lange et al. 2017)
              crohn's disease, (de Lange et al. 2017)
              juvenile arthtiris, (López-Isac et al 2020)
              primary sclerosing cholangitis (PSC), 
              type 1 diabetes, (Ji et al 2016)
              rehumatoid artritis, (the one the okada sent me via mail, all info in the excell) (Okada et al. 2014)
              asthma, (Han et al 2020)
              atopicdermatitis. (Sliz et al 2021)
              systemic lupus erythematosus (Bentham et al 2015)
The code (GitLab) for this GWAS is in brach uc-cd-psc-ra-sle-t1d-jia-asthm-derma.
All info about the summary statistics is in the excell.file in our Teams channel (AID_sharing) and soon uploaded here. 
          
          
e) 2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/ is the SAME AS gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/, but we change the sample size of asthma as it was a mistake in GWAS catalog. Moreover, the MAF from the GWAS was used when available. 
it was a 3 factor GWAS with the following diseases: 
              ulcerative colitis, (de Lange et al. 2017)
              crohn's disease, (de Lange et al. 2017)
              juvenile arthtiris, (López-Isac et al 2020)
              primary sclerosing cholangitis (PSC), 
              type 1 diabetes, (Ji et al 2016)
              rehumatoid artritis, (the one the okada sent me via mail, all info in the excell) (Okada et al. 2014)
              asthma, (Han et al 2020)
              atopicdermatitis. (Sliz et al 2021)
              systemic lupus erythematosus (Bentham et al 2015)
The code (GitLab) for this GWAS is in brach uc-cd-psc-ra-sle-t1d-jia-asthm-derma.
All info about the summary statistics is in the excell.file in our Teams channel (AID_sharing) and soon uploaded here. 
  
  
f) 3_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma_NO-MAF-INFO/ is the SAME AS 2_gwas_uc-cd-psc-ra-sle-t1d-jia-asthma-derma/, (so the sample size of asthma is the correcte one), but we removed the MAF and INFO column from asthma_han. This  because in the other two versions the MAF column and INFO columns in asthma was used to filter out SNPs in asthma_han, but in all the other GWAS the MAF column was not there (or was not interpreted by genomic sem) and therefore the MAF filtering was applyed only for asthma. Here all of the traits were thus prepared equally using the filtering that the package does with the reference file. 
it was a 3 factor GWAS with the following diseases: 
              ulcerative colitis, (de Lange et al. 2017)
              crohn's disease, (de Lange et al. 2017)
              juvenile arthtiris, (López-Isac et al 2020)
              primary sclerosing cholangitis (PSC), 
              type 1 diabetes, (Ji et al 2016)
              rehumatoid artritis, (the one the okada sent me via mail, all info in the excell) (Okada et al. 2014)
              asthma, (Han et al 2020)
              atopicdermatitis. (Sliz et al 2021)
              systemic lupus erythematosus (Bentham et al 2015)
The code (GitLab) for this GWAS is in brach uc-cd-psc-ra-sle-t1d-jia-asthm-derma.
All info about the summary statistics is in the excell.file in our Teams channel (AID_sharing) and soon uploaded here. 