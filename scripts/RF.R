library(tidymodels)
library(RSQLite)
library(future)

plan(multisession)
tidymodels_prefer()

splits <- read.table("data/assignment_sets.txt", header = T)
marker <- read.table("data/Marker_matrix.csv", header = T, sep = ";"); m_info <- "snp"

y <- read.table("data/hs_br.txt", header = T, sep = ";"); y_info <- "br"
# y <- read.table("data/hs_fus1.txt",      sep = ";", header = T);  y_info = "fus1"
# y <- read.table("data/hs_septoria.txt",  sep = ";", header = T);  y_info = "septoria"
# y <- read.table("data/hs_md.txt",        sep = ";", header = T);  y_info = "md"
# y <- read.table("data/hs_yr.txt",        sep = ";", header = T);  y_info = "yr"

data <- cbind(y = y$yy, marker)

rfile  <- paste(y_info, m_info)
dbfile <- "RF.sqldb"
add.df <- data.frame(method="rf", data.x=m_info, data.y=0, run=0, trait=y_info, cor=0, time=0)

i=1
for(i in 1:200){
    idx <- splits[,i] == "T"
    train <- data[idx,]
    test  <- data[!idx,]

    train.cv <- vfold_cv(train,5)

    rf_model <- rand_forest(
        mtry = tune(),
        trees = tune(),
        min_n = tune(),
        mode = "regression") %>% set_engine(engine = "ranger",num.threads = 10)

    t0 <- Sys.time()
    rf_recep <- recipe(train) %>%
        update_role(everything(), new_role = "predictor") %>% 
        update_role(y, new_role = "outcome")

    rf_wflow <- workflow(preprocessor = rf_recep, spec = rf_model)

    # Hyperparameter space
    rf_set <- extract_parameter_set_dials(rf_wflow)
    rf_set <-  rf_set %>% 
        update(mtry  = mtry( c(round(ncol(train)*0.01), round(ncol(train)/3)))) %>% 
            update(min_n = min_n(c( 1,20) ) ) %>% 
                update(trees = trees(c(200,1000)) )

    # Entropy grid generates 10 hyperparameter combinations for training
    entropy.grid <- grid_max_entropy(rf_set, size = 10)
    rf_initial <-  rf_wflow %>%
     tune_grid(resamples = train.cv,
               grid = entropy.grid,
               metrics = metric_set(rmse))
    
    # Bayesian search samples 10 more hyperparameter combinations to train
    bayes.res <- tune_bayes(
      object = rf_wflow,
      resamples = train.cv,
      param_info = rf_set,
      initial = rf_initial,
      iter = 10,
      metrics = metric_set(rmse),
      control = control_bayes(verbose = TRUE)
    )
    
    final_model <-
        finalize_workflow(rf_wflow, select_best(bayes.res, metric = "rmse")) %>%
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
    gridResults$algo <- "RF"

    # Write everything in database: correlation, time required, predicted and original values
    # Also save results from random grid search
    conn <- dbConnect(RSQLite::SQLite(), dbfile)
    sqliteSetBusyHandler(conn, 30000)
    dbWriteTable(conn, rfile, add.df, append=TRUE)
    dbWriteTable(conn, paste0(rfile,"_pred"), pred.t, append=TRUE)
    dbWriteTable(conn, paste0(rfile,"_orig"), orig.t, append=TRUE)
    dbWriteTable(conn, paste0(rfile,"_grid"), gridResults, append=TRUE)
    dbDisconnect(conn)
}