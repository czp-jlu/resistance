#library(parsnip)
library(tidymodels)
library(bonsai)
library(RSQLite)

library(future)
plan(multisession)

splits <- read.table("data/assignment_sets.txt", header = T)
marker <- read.table("data/Marker_matrix.csv", header = T, sep = ";"); m_info <- "snp"

y <- read.table("data/hs_br.txt",        sep = ";", header = T); y_info <- "br"
# y <- read.table("data/hs_fus1.txt",      sep = ";", header = T);  y_info = "fus1"
# y <- read.table("data/hs_septoria.txt",  sep = ";", header = T);  y_info = "septoria"
# y <- read.table("data/hs_md.txt",        sep = ";", header = T);  y_info = "md"
# y <- read.table("data/hs_yr.txt",        sep = ";", header = T);  y_info = "yr"

data <- cbind(y = y$yy, marker)

rfile  <- paste(y_info, m_info,"GB")
dbfile <- "GB.sqldb"
add.df <- data.frame(method="gbm_c", data.x=m_info, data.y=0, run=0, trait=y_info, cor=0, time=0)

tidymodels_prefer()


data2 <- data
threshold <- quantile(data2$y, .1)
hist(data2$y)

data2$y <- ifelse(data2$y <= threshold, "Top","Flop")
data2$y <- as.factor(data2$y)
table(data2$y)

for(i in 1:200){
  
  idx <- splits[,i] == "T"
  train <- data2[idx,]
  test  <- data2[!idx,]
  
  train.cv <- vfold_cv(train,5,strata = "y")
  
  t0 <- Sys.time()
  gb_model <- boost_tree(
    trees = tune(),
    min_n = 1,
    mtry = tune(),
    learn_rate = .01,
    mode = "classification"
  ) %>% set_engine("lightgbm", counts = F, num_threads = 10)
  
  gb_recep <- recipe(train) %>%
    update_role(everything(colnames(train)[2:16668]), new_role = "predictor") %>% 
    update_role(y, new_role = "outcome")  
  
  gb_wflow <- workflow(preprocessor = gb_recep, spec = gb_model)
  
  gb_set <- extract_parameter_set_dials(gb_wflow)
  
  # Hyperparameter space
  gb_set <- gb_set %>%
    update(mtry = mtry_prop(c(0.01,0.8)) ) %>%
        update(trees = trees(c(500,2000)) )
  
  # Entropy grid generates 10 hyperparameter combinations for training
  entropy.grid <- grid_max_entropy(gb_set, size = 10)
  
  gb_initial <-  gb_wflow %>%
    tune_grid(resamples = train.cv,
              grid = entropy.grid,
              metrics = metric_set(kap))
  
  # Bayesian search samples 10 more hyperparameter combinations to train
  bayes.res <- tune_bayes(
    object = gb_wflow,
    resamples = train.cv,
    param_info = gb_set,
    initial = gb_initial,
    iter = 10,
    metrics = metric_set(kap),
    control = control_bayes(verbose = TRUE)
  )
  
  # Not sure if automatically considers the entropy grid too when
  # selecting the model with the highest kappa,
  # therefore we make sure it does
  eg1 <- show_best(gb_initial,1,metric = "kap")
  ba1 <- show_best(bayes.res, 1,metric = "kap")

  if(eg1$mean < ba1$mean){
      best <- select_best(bayes.res, "kap")} else {
          best <- select_best(gb_initial, "kap")}

  final_model <- finalize_workflow(gb_wflow, best) %>%
    fit(train)

  pr <- manual_final %>% predict(test)

  t1 <- Sys.time()
  
  #  This will be appended to the database, contains all infos of the run
  add.df$cor <- kap_vec(pr$.pred_class, test$y)
  add.df$run <- i
  add.df$time <- as.numeric(t1-t0)
  
  
  # Prepare data for saving
  gridResults <- as.data.frame(do.call(rbind, bayes.res$.metrics))
  pred.t <- as.data.frame(t(pr))
  orig.t <- as.data.frame(t(test$y))
  gridResults$run <- i
  gridResults$m_type <- m_info
  gridResults$algo <- "GB"
  
  # Write everything in database: correlation, time required, predicted and original values
  # Also save results from random grid search
  conn <- dbConnect(RSQLite::SQLite(), dbfile)
  sqliteSetBusyHandler(conn, 30000)
  dbWriteTable(conn, rfile, add.df, append=TRUE)
  dbWriteTable(conn, paste0(rfile," pred"), pred.t, append=TRUE)
  dbWriteTable(conn, paste0(rfile," orig"), orig.t, append=TRUE)
  dbWriteTable(conn, paste0(rfile," grid"), gridResults, append=TRUE)
  dbDisconnect(conn)
}
