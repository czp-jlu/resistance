#######################################
# HaploSelekt: Cross-validation study #
# Carola Zenke-Philippi               #
#######################################
# Prediction of resistance scores     #
# from marker data with BGLR          #
# (Bayesian ordinal regression)       #
#######################################

library("BGLR")
library("RSQLite")
library("SelectionTools")

# data base files
dbfile <- "sqldb/hs_bglr.sqldb" # database
rfile.cor <- "results.cor" # file for the results in the database (correlation)
rr <- data.frame()        # prepare data frame for the results

# other information
NRUNS <- 1:200 # 200 runs

# assignment of genotypes to training and validation set
sets <- read.table("data/assignment_sets.txt", header=T) # file contains information for 1000 runs
sets[1:5,1:5]

# traits
traits <- c("br", "yr", "fus1", "md", "septoria")

# prepare marker data
st.read.marker.data("data/haploselekt-beagle-out.mpo", format="m",
                    data.set="default") # data need to be imputed with BEAGLE first
ZZ <- gs.build.Z(out.filename="ZZ", auxfiles=T, data.set="default")

# parameters for the model
ETA <- list( list( X=ZZ, model="BRR" ) )

# cross validation
for (RUN in NRUNS) {
  
  set <- as.vector(sets[,RUN])
  
  for (TRT in traits) {
    
    # read phenotypic data
    fname <- paste0("data/hs_", TRT, ".txt")
    yy <- read.table(fname, header=T, sep=";")
    
    yy.ts <- yy$yy
    yy.ts[set=="P"] <- NA
    yy.ps <- yy$yy
    yy.ps <- yy.ps[set=="P"]
    
    t0 <- Sys.time()
    
    # fit model 
    AA <- BGLR(y=yy.ts, 
               response_type="ordinal", 
               ETA=ETA, 
               nIter = 1e4, 
               burnIn = 1e3, 
               verbose=F,
               saveAt = "bglr/")
    
    # predicted values
    yhat <- AA$yHat[set=="P"]
    
    t1 <- Sys.time()
    
    tt <- as.numeric(t1-t0)
    
    # observed values
    yy.result <- data.frame(matrix(yhat, byrow=T, nrow=1))
    rfile.yy <- paste0("yhat_bglr_snp_0_", TRT )
    
    conn <- dbConnect(RSQLite::SQLite(), dbfile)
    dbWriteTable(conn, rfile.yy, yy.result, append=TRUE)
    dbDisconnect(conn)
    
    # prediction accuracy
    cc <- cor(yhat, yy.ps)

    # save results
    rr.cv <- data.frame(method = "bglr",
                        data.x = "snp",
                        data.y = 0,
                        run = RUN,
                        trait = TRT,
                        cor = cc,
                        time = tt)
    
    # Save simulation results in data base
    conn <- dbConnect(RSQLite::SQLite(), dbfile)
    dbWriteTable(conn, rfile.cor, rr.cv, append=TRUE)
    dbDisconnect(conn)
    
  }
}
