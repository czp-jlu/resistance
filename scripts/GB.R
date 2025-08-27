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
add.df <- data.frame(method="gbm", data.x=m_info, data.y=0, run=0, trait=y_info, cor=0, time=0)

tidymodels_prefer()

i=1
for(i in 1:200){
  idx <- splits[,i] == "T"
  train <- data[idx,]
  test  <- data[!idx,]
  
  train.cv <- vfold_cv(train,5)

  t0 <- Sys.time()
  gb_model <- boost_tree(
      trees = tune(),
      min_n = tune(),
      mtry = tune(),
      learn_rate = .001,
      mode = "regression"
  ) %>% set_engine("lightgbm", counts = F, num_threads = 10)
  
  gb_recep <- recipe(train) %>%
    update_role(everything(), new_role = "predictor") %>% 
    update_role(y, new_role = "outcome")
  
  gb_wflow <- workflow(preprocessor = gb_recep, spec = gb_model)
  
  gb_set <- extract_parameter_set_dials(gb_wflow)

  # Hyperparameter space
  gb_set <- gb_set %>%
    update(mtry = mtry_prop(c(0.01,0.2)) ) %>%
    update(trees = trees(c(50,500)) )

  # Entropy grid generates 10 hyperparameter combinations for training
  entropy.grid <- grid_max_entropy(gb_set, size = 10)
  gb_initial <-  gb_wflow %>%
     tune_grid(resamples = train.cv,
               grid = entropy.grid,
               metrics = metric_set(rmse))
  
  # Bayesian search samples 10 more hyperparameter combinations to train
  bayes.res <- tune_bayes(
    object = gb_wflow,
    resamples = train.cv,
    param_info = gb_set,
    initial = gb_initial,
    iter = 10,
    metrics = metric_set(rmse),
    control = control_bayes(verbose = TRUE)
  )

  final_model <-
    finalize_workflow(gb_wflow, select_best(bayes.res, metric = "rmse")) %>%
    fit(train)

  pr <- final_model %>% predict(test)
  t1 <- Sys.time()
  
  #  This will be appended to the database, contains all infos of the run
  add.df$cor <- cor(pr$.pred, test$y)
  add.df$run <- i
  add.df$time <- as.numeric(difftime(t1, t0, units = 'mins'))


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
  sqliteSetBusyHandler(conn, 9000)
  dbWriteTable(conn, rfile, add.df, append=TRUE)
  dbWriteTable(conn, paste0(rfile," pred"), pred.t, append=TRUE)
  dbWriteTable(conn, paste0(rfile," orig"), orig.t, append=TRUE)
  dbWriteTable(conn, paste0(rfile," grid"), gridResults, append=TRUE)
  dbDisconnect(conn)
}