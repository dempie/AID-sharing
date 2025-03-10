# AID_sharing project

**Repository for the paper**: [https://www.nature.com/articles/s41467-023-38389-6](https://www.nature.com/articles/s41467-023-38389-6)

**Contributors**: Pietro Demela, Nicola Pirastu, Blagoje Soskic

The final version of the code (February 2023) that was used to generate the figures and analyses for the Nat Com paper is in the `revision_1` branch, in the `factor_gwas` folder.

---

## Project Aim

The aim of the project is to use structural equation modelling to model the polygenic structure of autoimmune diseases. The package used is **GenomicSEM**, which can be found here: [https://github.com/GenomicSEM/GenomicSEM](https://github.com/GenomicSEM/GenomicSEM).

---

## Home Directory Structure

This is the home directory of the project. Below is a description of what is contained in each folder:

### a) `ldscores/`
- A folder downloaded as described by the GitHub of **GenomicSEM**, that contains the LD scores from HapMap3.
- **DO NOT MODIFY ANYTHING**.

### b) `outputs/`
- All outputs will be inside this folder in subfolders according to the version.

### c) `Prevalences/`
- This folder contains CSV files generated by reading the papers of the GWAS.
- Each CSV contains the sample prevalence for each cohort in the GWAS. 
- It's necessary for calculating the sample effective sample size. Details are in the GitHub of the package.

### d) `renv/`
- The folder generated by the **renv** package to keep track of the package versions.

### e) `Scripts/`
- Contains the code written by the contributors, where the analysis is conducted.

### f) `SNP/`
- Contains a reference file of the rsID needed by the package.
- This file is also used to add rsID to summary stats.
- Information on how to use it is available in the GitHub of **GenomicSEM**.

### g) `Summary_Stats/`
- Contains the raw summary stats downloaded from the papers, mostly from **GWAS catalog**.


   
