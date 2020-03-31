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

# lasso model
bc.fit.lasso <- glmnet(x.bc.dt, y.bc.dt, subset = train.idx, family="binomial")  # same as setting alpha=1

# ridge regression
bc.fit.ridge <- glmnet(x.bc.dt, y.bc.dt, alpha=0, subset = train.idx, family="binomial")

# plotting the trajectories of the coefficients for various lambda λ
par(mfrow=c(1,2), mar=c(4,4,5,2))
plot(bc.fit.lasso, main="Lasso trajectories") 
plot(bc.fit.ridge, main="Ridge trajectories")

# for each model learn by cross-validation the penalty parameter λ that maximizes the AUC

# lasso model
bc.fit.cv.lasso <- cv.glmnet(x.bc.dt, y.bc.dt, subset = train.idx, family="binomial", type.measure="auc")

# ridge regression
bc.fit.cv.ridge <- cv.glmnet(x.bc.dt, y.bc.dt, subset = train.idx, family="binomial", type.measure="auc", alpha=0)

# plotting the cross-validation curve

par(mfrow=c(1,2), mar=c(4,4,5,2)) 
plot(bc.fit.cv.lasso, main="Lasso") 
plot(bc.fit.cv.ridge, main="Ridge")

# penalty parameters λ that maximizes the AUC
bc.fit.cv.lasso$lambda.min # this is admittedly counterintuitive
bc.fit.cv.ridge$lambda.min # this is admittedly counterintuitive
