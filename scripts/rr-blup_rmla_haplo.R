#######################################
# HaploSelekt: Cross-validation study #
# Carola Zenke-Philippi               #
#######################################
# Predictions with RR-BLUP/RMLA       #
# with haplotype blocks (V28)         #
# for untransformed and logit data    #
#######################################

library ("SelectionTools")
library("sqldf")

NRUNS <- 1:200 # 200 runs

# data base
dbfile <- "sqldb/hs_rrblup_haplo.sqldb" # database
rfile.cor <- "results.cor" # file for the results in the database (correlation)
rr <- data.frame()        # prepare data frame for the results

# set input and output folders
st.input.dir  <- "input"
st.output.dir <- "output"

# create output folder
dir.create(st.output.dir)

# traits
traits <- c("br", "yr", "fus1", "md", "septoria")

###################################
# Copy files for the versions to the
# input directory:

if (!("input" %in% list.files())) {
  
  dir.create(st.input.dir)
  
}

files.list <- c("haploselekt-map.txt", "haploselekt-beagle-out.mpo")
# marker data must be created with BEAGLE first
file.copy(from=paste0("data/", files.list), to=paste0(st.input.dir, "/", files.list), overwrite=T)

# files for the phenotypic data
files.list <- c(paste0("hs_", traits, ".txt"), paste0("hs_", traits, "_logit.txt"))
file.copy(from=paste0("data/", files.list), to=paste0(st.input.dir, "/", files.list), overwrite=T)

###################################
# Different assignment of lines to the
# training and validation set:

sets <- read.table("data/assignment_sets.txt", header=T) # sets for 1000 runs

###################################
# Parameters that are the same for all
# CV runs:

rr <- data.frame()        # prepare data frame for the results
nn <- 361                 # number of genotypes in the preselected data set

gs.set.num.threads(4)
st.set.info.level (-2)

# file names
fname1 <- paste0(st.input.dir, "/hs-haplo-") # marker data
fname2 <- "hs-haplo-" # marker data

# estimation methods
estimation <- c("rr-blup", "rmla")

# tranformation
transformation <- c(0, 1) # 0: untransformed y values, 1: logit-transformed y values

###################################
# marker data

marker <- read.table("input/hs-hblocks.txt", header=T)
# tutorial on creation of haplotype blocks at
# https://population-genetics.uni-giessen.de/~software/
# or see supplementary files on Difabachew et al. 2023

###################################
# ACTUAL CROSS VALIDATION
###################################

for (RUN in NRUNS) {
  
  # validation and training set
  no.v.set <- c(1:nn)[sets[,RUN]=="P"]
  no.t.set <- c(1:nn)[sets[,RUN]=="T"]
  
  v.set <- marker[,no.v.set] # markers only for individuals in the validation set
  t.set <- marker[,no.t.set] # markers only for individuals in the training set
    
  write.table(v.set, paste0(fname1, "-v.mpo"), quote=F)
  write.table(t.set, paste0(fname1, "-t.mpo"), quote=F)
    
  st.read.marker.data (paste0(fname2, "-v.mpo"), # load marker data for prediction set
                         format="m",
                         data.set = "vv")
  st.read.marker.data (paste0(fname2, "-t.mpo"), # load marker data for training set
                         format="m",
                         data.set = "tt")
    
  # Preprocessing of data
  st.restrict.marker.data(NoAll.MAX=2, data.set="tt") # Maximum number of alleles
  st.restrict.marker.data(MaMis.MAX=0.2, data.set="tt") # Max missing at a marker
  st.restrict.marker.data(ExHet.MIN=0.05, data.set="tt") # Minimum gene diversity      
    
  for (TRT in traits){
    
    fnames <- c(paste0("hs_", TRT, ".txt"), paste0("hs_", TRT, "_logit.txt"))
    
    for (TRS in transformation) {
      
      fname3 <- fnames[TRS + 1]
  
      # Read in phenotypic data
      st.read.performance.data(fname3, data.set="tt") # phenotypic data for TS
      st.read.performance.data(fname3, data.set="vv") # phenotypic data for VS
        
        for (EST in estimation){
          
          t0 <- Sys.time()
          
          if (EST == "rr-blup"){
            gs.esteff.rr ( method="BLUP", data.set="tt" )
          } else {
            gs.esteff.rr ( method="RMLA", data.set="tt" )
          } # is/else EST
          
          phe <- gs.predict.genotypes ( training.set   = "tt", # prediction of phenotypes in the VS
                                        prediction.set = "vv")
          
          t1 <- Sys.time()
          
          tt <- as.numeric(t1-t0)

          if ( TRS == 1 ) {

              phe$yhat <- (exp(phe$yhat)/(1+exp(phe$yhat)))*10

          }
  
          # observed values
          yy.result <- data.frame(matrix(phe$yhat, byrow=T, nrow=1))
          rfile.yy <- paste0("yhat_", EST, "_hb_", TRS, "_", TRT )
          
          conn <- dbConnect(RSQLite::SQLite(), dbfile)
          dbWriteTable(conn, rfile.yy, yy.result, append=TRUE)
          dbDisconnect(conn)
          
          # prediction accuracy
          cc <- cor(phe$y, phe$yhat)
          
          # save results
          rr.cv <- data.frame(method = EST,
                              data.x = "hb",
                              data.y = TRS,
                              run = RUN,
                              trait = TRT,
                              cor = cc,
                              time = tt)
          
          # Save simulation results in data base
          conn <- dbConnect(RSQLite::SQLite(), dbfile)
          dbWriteTable(conn, rfile.cor, rr.cv, append=TRUE)
          dbDisconnect(conn)
          
          } # EST
      } # for TRS
    } # for TRT
  }# for RUN
