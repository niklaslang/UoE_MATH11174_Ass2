### (a) ###

library(data.table)
bc.dt <- fread("wdbc2.csv", stringsAsFactors = TRUE)

# Using package caret, create a data partition so that the training set 
# contains 70% of the observations (set the random seed to 1 beforehand)

library(caret)

set.seed(1)
train.idx <- createDataPartition(bc.dt$diagnosis, p=0.7)$Resample1

# Fitting both a ridge regression model and a lasso model on the training set 
# to diagnose the type of tumour from the 30 biomarkers

library(glmnet)

# function to transform a dataframe to a matrix as expected by the glmnet package

prepare.glmnet <- function(data, formula=~ .){
  
  ## create the design matrix to deal correctly with factor variables, ## without losing rows containing NAs
  old.opts <- options(na.action='na.pass')
  x <- model.matrix(formula, data)
  options(old.opts)
  
  ## remove the intercept column, as glmnet will add one by default 
  x <- x[, -match("(Intercept)", colnames(x))]
  return(x)
}

x.bc.dt <- prepare.glmnet(bc.dt[,!"id"], ~ . - diagnosis) # exclude the outcome
y.bc.dt <- bc.dt$diagnosis # store the outcome separately

# ridge regression
bc.ridge <- glmnet(x.bc.dt, y.bc.dt, alpha=0, subset = train.idx, family="binomial")
# lasso model
bc.lasso <- glmnet(x.bc.dt, y.bc.dt, subset = train.idx, family="binomial")  # same as setting alpha=1

# plotting the trajectories of the coefficients for various lambda λ
par(mfrow=c(1,2), mar=c(4,4,5,2))
plot(bc.ridge, main="Ridge trajectories")
plot(bc.lasso, main="Lasso trajectories") 

# for each model learn by cross-validation the penalty parameter λ that maximizes the AUC
# ridge regression
bc.cv.ridge <- cv.glmnet(x.bc.dt, y.bc.dt, subset = train.idx, family="binomial", type.measure="auc", alpha=0)
# lasso model
bc.cv.lasso <- cv.glmnet(x.bc.dt, y.bc.dt, subset = train.idx, family="binomial", type.measure="auc")

# plotting the cross-validation curve

par(mfrow=c(1,2), mar=c(4,4,5,2)) 
plot(bc.cv.ridge, main="Ridge")
plot(bc.cv.lasso, main="Lasso") 


# penalty parameters λ that maximizes the AUC
bc.cv.ridge$lambda.min # this is admittedly counterintuitive
bc.cv.lasso$lambda.min # this is admittedly counterintuitive

### (b) ###

# for both models fitted in (a) extract the AUCs corresponding to the optimal λ 

# There are multiple ways for retrieve the AUC obtained with optimal λ:

# option 1: optimal λ maximzes the AUC, so we just extract the maximal AUC:
# ridge
max(bc.cv.ridge$cvm)
# lasso
max(bc.cv.lasso$cvm)

# option 2: extract the AUC from the field cvm corresponding to the optimal λ (lamda.min) from the field lambda

# ridge
bc.cv.ridge$cvm[which(bc.cv.ridge$lambda == bc.cv.ridge$lambda.min)]
# lasso
bc.cv.lasso$cvm[which(bc.cv.lasso$lambda == bc.cv.lasso$lambda.min)]

# and to the λ such that the AUC is within 1 standard error of the maximum
# ridge
bc.cv.ridge$cvm[which(bc.cv.ridge$lambda <= bc.cv.ridge$lambda.1se)]
# lasso
bc.cv.lasso$cvm[which(bc.cv.lasso$lambda <= bc.cv.lasso$lambda.1se)]

