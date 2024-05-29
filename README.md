This repository contains the genotypic and phenotypic data as well as code for the manuscript "Machine learning for prediction of resistance scores in wheat".

#### Data:

data/assignment_sets.txt: Assignment of genotypes to training and validation sets. <br>
data/haploselekt-map.txt: Genetic map. <br>
data/haploselekt-marker.mpo: Unfiltered marker matrix. <br>
data/hs_br.txt: Resistance scores for brown rust (P. triticina). <br>
data/hs_fus1.txt: Resistance scores for Fusarium (F. graminearum). <br>
data/hs_md.txt: Resistance scores for mildew (B. graminis). <br>
data/hs_septoria.txt: Resistance scores for Septoria (S. graminis). <br>
data/hs_yr.txt: Resistance scores for yellow rust (P. striiformis). <br>
data/Marker_matrix.csv: Marker matrix used for the machine learning methods. <br>
input/hs-blocks.txt: File with haplotype blocks.

#### Scripts:

Autoencoder.py: Python script to create the autoencoder data. <br>
beagle.R: R script to impute the missing marker data with BEAGLE. <br>
bglr.R: R script for method BGLR. <br>
GB.R: R script for method GB (regression). <br>
GB_classif.R: R script for method GB (classification). <br>
RF.R: R script for method RF. <br>
rr-blup_rmla.R: R script for methods RR-BLUP and RMLA with single SNPs. <br>
rr-blup_rmla_haplo.R: R script for methods RR-BLUP and RMLA with haplotype blocks. <br>
SVM_master.R: R script to run multiple instances of SVM_worker.R or SVM_worker_classif.R at the same time. <br>
SVM_worker.R: R script for SVM (regression). Requires arguments passed through SVM_master.R. <br>
SVM_worker_classif.R: R script for SVM (classification). Requires arguments passed through SVM_master.R. <br>
transformation.R: R script for the logit transformation of resistance scores.

For the incremental feature selection, please refer to [Heinreich et al. 2023](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-023-00853-8) and their github repository [here](https://github.com/FelixHeinrich/GP_with_IFS/)
