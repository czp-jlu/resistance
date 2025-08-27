library(parallel)

runner <- function(args){
  system(paste('Rscript --max-ppsize=500000 $EDIT_FULL_PATH_TO_YOUR_SCRIPT_HERE$/SVM_worker.R',
               args$i, args$name, args$trait, getwd()))
  }

# runner <- function(args){
#   system(paste('Rscript --max-ppsize=500000 $EDIT_FULL_PATH_TO_YOUR_SCRIPT_HERE$/SVM_worker_classif.R',
#                args$i, args$name, args$trait, getwd()))
# }

merger <- function(path="SVM_res"){
  all.files <- list.files(path)
  all.files.cor <- all.files[grep("cor", all.files)]
  all.files.cor <- paste(path,all.files.cor,sep="/")
  
  all.files.pred <- all.files[grep("pred", all.files)]
  all.files.pred <- paste(path,all.files.pred,sep="/")
  
  cors <- lapply(all.files.cor, read.table, sep=";")
  cors <- do.call(rbind, cors)
  file.remove(all.files.cor)
  
  pred <- lapply(all.files.pred, read.table, sep=";")
  pred <- do.call(rbind, pred)
  file.remove(all.files.pred)
  
  write.table(cors,"cors_SVM")
  write.table(pred,"pred_SVM")
  
  unlink(path, recursive = T,force = T)
}

# 50 runs in parallel
num_cores <- 50
cl <- makeCluster(num_cores)
i.s <- seq(0,200,num_cores)

data.info <-
  data.frame(
    name  = c(
      "hs_br.txt",
      "hs_md.txt",
      "hs_yr.txt",
      "hs_sept.txt",
      "hs_fus1.txt"
    ),
    trait = c("br", "md", "yr", "septoria", "fus1"),
    time = 0
  )

for(k in 1:nrow(data.info)){
  name  <- data.info$name[k]
  trait <- data.info$trait[k]
  
  t0 <- Sys.time()
  for(i in 1:(length(i.s) - 1)) {
    start <- i.s[i] + 1
    end <- i.s[i + 1]
    
    df <- data.frame(i = start:end, name = name, trait = trait)
    args <- apply(df, 1, as.list)
    
    parLapply(cl, args , runner)
    
    gc()
  }
  t1 <- Sys.time()
  
  data.info$time[k] <- difftime(t1,t0, units="min")/200
  
}

# Close the cluster
stopCluster(cl)
write.table(data.info, "Results_SVR.csv")

merger()