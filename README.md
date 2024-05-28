This repository contains the genotypic and phenotypic data as well as code for the manuscript "Machine learning for prediction of resistance scores in wheat".

data/assignment_sets.txt: Assignment of genotypes to training and validation sets.
data/haploselekt-map.txt: Genetic map.
data/haploselekt-marker.mpo: Unfiltered marker matrix.
data/hs_br.txt: Resistance scores for brown rust (P. triticina).
data/hs_fus1.txt: Resistance scores for Fusarium (F. graminearum).
data/hs_md.txt: Resistance scores for mildew (B. graminis).
data/hs_septoria.txt: Resistance scores for Septoria (S. graminis).
data/hs_yr.txt: Resistance scores for yellow rust (P. striiformis).
data/Marker_matrix.csv: Marker matrix used for the machine learning methods.
input/hs-blocks.txt: File with haplotype blocks.

Scripts:
Autoencoder.py: Python script to create the autoencoder data.
beagle.R: R script to impute the missing marker data with BEAGLE.
bglr.R: R script for method BGLR.
GB.R: R script for method GB (regression).
GB_classif.R: R script for method GB (classification).
RF.R: R script for method RF.
rr-blup_rmla.R: R script for methods RR-BLUP and RMLA with single SNPs.
rr-blup_rmla_haplo.R: R script for methods RR-BLUP and RMLA with haplotype blocks.
SVM_master.R:
SVM_worker.R:
SVM_worker_classif.R:
transformation.R: R script for the logit transformation of resistance scores.
