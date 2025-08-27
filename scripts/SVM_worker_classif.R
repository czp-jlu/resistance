# Set options("expressions"=5e5)
# This is only relevant if your dataset exceeds a certain number of columns
# Refer to this stack overflow question for more information and possible 
# new updates
# https://stackoverflow.com/questions/78372245/tidymodels-workaround-for-kernlab-protection-stack-overflow-error

options("expressions"=5e5)

args <- commandArgs(trailingOnly = TRUE)

setwd(args[4])

suppressMessages(library(tidymodels))

tidymodels_prefer()

i=as.integer(args[1])

cat("Initialising replication ",i," ....")

splits <- read.table("data/assignment_sets.txt", header = T)
marker <- read.table("data/Marker_matrix.csv", header = T, sep = ";"); m_info <- "snp"

y <- read.table(args[2], header = T, sep = ";"); y_info <- args[3]

cat(" successful \n")

data <- cbind(y = y$yy, marker)
rfile  <- paste(y_info, m_info)

add.df <- data.frame(method="svm", data.x=m_info, data.y=0, run=0, trait=y_info, cor=0, time=0)

### Training sequence starts

idx <- splits[,i] == "T"

train <- data[idx,]
test  <- data[!idx,]

train.cv <- vfold_cv(train,5)

t0 <- Sys.time()

# requires R package 'kernlab' to be installed
svm_model <- svm_linear(
  cost = tune(),
  margin = tune(),
  engine = "kernlab",
  mode = "classification"
)

svm_recep <- recipe(x = train) %>%
update_role(everything(), new_role = "predictor") %>% 
update_role(y, new_role = "outcome")

svm_wflow <- workflow(preprocessor = svm_recep, spec = svm_model)

svm_set <- extract_parameter_set_dials(svm_wflow)

# Entropy grid generates 10 hyperparameter combinations for training
entropy.grid <- grid_max_entropy(svm_set, size = 10)
svm_initial <-  svm_wflow %>%
  tune_grid(resamples = train.cv,
            grid = entropy.grid,
            metrics = metric_set(kap))

# Bayesian search samples 10 more hyperparameter combinations to train
bayes.res <- tune_bayes(
  object = svm_wflow,
  resamples = train.cv,
  param_info = svm_set,
  initial = svm_initial,
  iter = 10,
  metrics = metric_set(kap),
  control = control_bayes(verbose = F)
)


final_model <- finalize_workflow(svm_wflow, select_best(bayes.res, metric = "kap")) %>%
  fit(train)

pr <- final_model %>% predict(test)
t1 <- Sys.time()

add.df$cor <- cor(pr$.pred, test$y)
add.df$run <- i
add.df$time1 <- t0
add.df$time2 <- t1

# Prepare data for saving
#gridResults <- as.data.frame(do.call(rbind, bayes.res$.metrics))
pred.t <- as.data.frame(t(pr))
orig.t <- as.data.frame(t(test$y))
#gridResults$run <- i
#gridResults$m_type <- m_info
#gridResults$algo <- "SVM"

# Databases may cause errors here due to seemingly non-functional busy handlers
# Instead, we create temporary files

if(!exists(paste0("SVM_res",m_info))){
  dir.create(paste0("SVM_res",m_info), showWarnings = F)
}

write.table(add.df, paste0("SVM_res/cor_",y_info,i), sep=";")
write.table(cbind(pred.t,i,m_info,y_info), paste0("SVM_res/pred_",y_info,i), sep=";")