#######################################
# HaploSelekt: Cross-validation study #
# Carola Zenke-Philippi               #
#######################################
# Logit-transformation of resistance  #
# scores                              #
#######################################

traits <- c("br", "yr", "fus1", "md", "septoria")

for (TRT in traits) {
  
  fname1 <- paste0("data/hs_", TRT, ".txt")
  fname2 <- paste0("data/hs_", TRT, "_logit.txt")
  
  dta <- read.table(fname1, header=T, sep=";")
  dta$yy <- dta$yy/10
  dta$yy <- log(dta$yy/(1-dta$yy))
  
  write.table(dta, fname2, quote=F, row.names=F, sep=";")
}
