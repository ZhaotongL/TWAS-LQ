# TWAS-LQ
This repository contains simulation and analysis codes for the paper "Accounting for nonlinear effects of gene expression identifies additional associated genes in transcriptome-wide association studies" (https://doi.org/10.1093/hmg/ddac015).

* ```TWASX2_bs_analysis_tede_lme_ukb.R``` is the script to perform TWAS-L, TWAS-Q, TWAS-LQ and their -LME versions, as well as TEDE test using UKB individual level data.
* ```TWASX2_bs_analysis_igap.R``` is the script to perform TWAS-L, TWAS-Q and TWAS-LQ using IGAP AD GWAS summary data.
* ```stage2.R``` is the function to perform TWAS stage 2 analysis based on disease/trait GWAS summary statistics.
* ```my_TEDE.R``` is the function to perform TEDE test extended to one or more terms of imputed GE.
* ```allele_qc.R``` is the function to align alleles between two studies.

