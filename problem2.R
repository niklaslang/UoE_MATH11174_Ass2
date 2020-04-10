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

### (c) ###

# Creating a data table for each model that reports 
# the choice of λ, the corresponding model size and their training AUC
# (3 significant digits for floating point values)

# ridge

bc.ridge.eval.dt <- data.table("λ" = signif(bc.cv.ridge$lambda,3),
                            "ModelSize" = bc.cv.ridge$nzero,
                            "AUC" = signif(bc.cv.ridge$cvm, 3))

# lasso 

bc.lasso.eval.dt <- data.table("λ" = signif(bc.cv.lasso$lambda,3),
                            "ModelSize" = bc.cv.lasso$nzero,
                            "AUC" = signif(bc.cv.lasso$cvm, 3))

# plotting the results

# ridge
par(mfrow=c(1,3), mar=c(4,4,5,2))
plot(bc.ridge.eval.dt$λ, bc.ridge.eval.dt$ModelSize, 
     main = "Ridge",
     xlab = 'λ',
     ylab = 'Number of non-zero parameters',
     col = "chartreuse3",
     cex = .5,
     ylim = c(0,30))
plot(bc.ridge.eval.dt$λ, bc.ridge.eval.dt$AUC, 
     main = "Ridge",
     xlab = 'λ',
     ylab = 'AUC',
     col = "aquamarine4",
     cex = .5,
     ylim = c(.95, 1))
plot(bc.ridge.eval.dt$ModelSize, bc.ridge.eval.dt$AUC, 
     main = "Ridge",
     xlab = 'Number of non-zero parameters',
     ylab = 'AUC', 
     col = "cadetblue3",
     cex = .5,
     xlim = c(0,30),
     ylim = c(.95, 1))

# lasso
par(mfrow=c(1,3), mar=c(4,4,5,2))
plot(bc.lasso.eval.dt$λ, bc.lasso.eval.dt$ModelSize, 
     main = "Lasso",
     xlab = 'λ',
     ylab = 'Number of non-zero parameters',
     col = "chartreuse3",
     cex = .5,
     ylim = c(0,30))
plot(bc.lasso.eval.dt$λ, bc.lasso.eval.dt$AUC, 
     main = "Lasso",
     xlab = 'λ',
     ylab = 'AUC',
     col = "aquamarine4",
     cex = .5,
     ylim = c(.95, 1))
plot(bc.lasso.eval.dt$ModelSize, bc.lasso.eval.dt$AUC, 
     main = "Lasso",
     xlab = 'Number of non-zero parameters',
     ylab = 'AUC', 
     col = "cadetblue3",
     cex = .5,
     ylim = c(.95, 1))

# commenting on results

### (d) ###

# Performing backwards elimination (we’ll later refer to this as model B) 
# on the same training set derived at point (a)

# standardize variables to obtain standardized coefficients
bc.sd1.dt <- copy(bc.dt)
covar.variables <- colnames(bc.dt[, -c("id", "diagnosis")]) # exclude non covariate columns
bc.sd1.dt[, (covar.variables) := lapply(.SD, function(x) x / sd(x, na.rm = TRUE)), .SDcols = covar.variables]

# load library for stepwise regression
library(MASS)

# fit logistic regression model with all the training split from (a)
bc.full.model <- glm(diagnosis ~ ., data=bc.sd1.dt[, !"id"], subset = train.idx, family="binomial") # all variables except the patient id

bc.sel.back <- stepAIC(bc.full.model, direction="back")

# Reporting the variables selected and their standardized regression coefficients 
# in decreasing order of the absolute value of their standardized regression coefficient

bc.sel.back.coefs <- data.table("variable" = rownames(coef(summary(bc.sel.back)))[-1],
                               "stdz.coef" = signif(coef(summary(bc.sel.back))[-1,1]))

#bc.el.back.coefs <- bc.el.back.coefs[variable != "(Intercept)"]
bc.sel.back.coefs <- bc.sel.back.coefs[order(-abs(bc.sel.back.coefs$stdz.coef))]
bc.sel.back.coefs

### (e) ###

# Repeating the same analysis of point (d) by using stepwise selection (model S)
# starting from the null model

bc.train.dt <- bc.sd1.dt[, !"id"][train.idx]

bc.null.model <- glm(diagnosis ~ 1, data=bc.sd1.dt[, !"id"], subset = train.idx, family="binomial")

bc.sel.forw <- step(bc.null.model, scope=list(upper=bc.full.model),direction="both")

# Reporting the variables selected and their standardized regression coefficients 
# in decreasing order of the absolute value of their standardized regression coefficient

bc.sel.forw.coefs <- data.table("variable" = rownames(coef(summary(bc.sel.forw)))[-1],
                                "stdz.coef" = signif(coef(summary(bc.sel.forw))[-1,1]))

bc.sel.forw.coefs <- bc.sel.forw.coefs[order(-abs(bc.sel.forw.coefs$stdz.coef))]
bc.sel.forw.coefs

# Did at any point in this procedure occur that a variable entered the model and was later on discarded? 
# If so which?

# Variables included in both models:
intersect(bc.sel.back.coefs$variable, bc.sel.forw.coefs$variable)

### (f) ###

# Comparing the goodness of fit of model `bc.sel.back`(backwards elimination) 
# and model `bc.sel.forw` (forward selection) using the AIC

# Given that the two models are not nested, we cannot use a likelihood ratio test, 
# but we can still use AIC and BIC to compare the models:

# AIC 
# takes into account model complexity: 
# the deviance of the model is penalised by twice the number of parameters estimated in the model
AIC(bc.sel.back)
AIC(bc.sel.forw)

# Model bc.sel.back has lower AIC, so it’s the better model

# BIC 
# BIC applies a stronger penalty for the use of additional predictors
BIC(bc.sel.back)
BIC(bc.sel.forw)

# in this case bc.sel.forw is better

## AUCs
#library(pROC)
#par(mfrow=c(1,1))
#roc(bc.dt$diagnosis[train.idx], bc.sel.back$fitted.values)
#roc(bc.dt$diagnosis[train.idx], bc.sel.forw$fitted.values)
#roc(bc.dt$diagnosis[train.idx], bc.sel.back$fitted.values, plot=TRUE, legacy.axes=TRUE, col="brown2", lwd = 1.4)
#roc(bc.dt$diagnosis[train.idx], bc.sel.forw$fitted.values, plot=TRUE, legacy.axes=TRUE, add=TRUE, col="deepskyblue3", lwd = 1.4)
#legend("bottomright", c("backwards elimination: AUC = 0.9916", "forward selection: AUC = 0.9891"), fill = c("brown2","deepskyblue3"))

### (g) ###

# Using model `bc.sel.back`` and model `bc.sel.forw` to predict observations in the training set, and from that compute
# the training AUC

bc.sel.back.pred <- predict(bc.sel.back, newdata=bc.dt[train.idx], type="response")
bc.sel.forw.pred <- predict(bc.sel.forw, newdata=bc.dt[train.idx], type="response")

par(mfrow=c(1,1))
roc(bc.dt$diagnosis[train.idx], bc.sel.back.pred)
roc(bc.dt$diagnosis[train.idx], bc.sel.forw.pred)
roc(bc.dt$diagnosis[train.idx], bc.sel.back$fitted.values, plot=TRUE, legacy.axes=TRUE, main = "Training AUCs", col="brown2", lwd = 3)
roc(bc.dt$diagnosis[train.idx], bc.sel.forw$fitted.values, plot=TRUE, legacy.axes=TRUE, add=TRUE, col="deepskyblue3", lwd = 3)
legend("bottomright", c("backwards elimination: AUC = 0.9916", "forward selection: AUC = 0.9891"), fill = c("brown2","deepskyblue3"))

### (h) ###

# Using the four models to predict the outcome for the observations in the test set 
# (use the lambda at 1 standard error for the penalised models)

# predict the outcome for the observations in the test data set
bc.sel.back.pred.obs <- predict(bc.sel.back, newdata=bc.dt[-train.idx,], type="response")
bc.sel.forw.pred.obs <- predict(bc.sel.forw, newdata=bc.dt[-train.idx,], type="response")
bc.ridge.pred.obs <- predict(bc.cv.ridge, newx = x.bc.dt[-train.idx,], type="response", s=bc.cv.ridge$lambda.1se)
bc.lasso.pred.obs <- predict(bc.cv.lasso, newx = x.bc.dt[-train.idx,], type="response", s=bc.cv.lasso$lambda.1se)

# Plotting the ROC curves of these models (on the same plot, using different colours) 
# and reporting their test AUCs. 

# compute test AUCs
roc(bc.dt$diagnosis[-train.idx], bc.sel.back.pred.obs)
roc(bc.dt$diagnosis[-train.idx], bc.sel.forw.pred.obs)
roc(bc.dt$diagnosis[-train.idx], as.vector(bc.ridge.pred.obs))
roc(bc.dt$diagnosis[-train.idx], as.vector(bc.lasso.pred.obs))

# plot test AUCs
par(mfrow=c(1,1))
roc(bc.dt$diagnosis[-train.idx], bc.sel.back.pred.obs, plot=TRUE, legacy.axes=TRUE, main = "Test AUCs",
    col="#FF9900FF", lwd = 3)
roc(bc.dt$diagnosis[-train.idx], bc.sel.forw.pred.obs, plot=TRUE, legacy.axes=TRUE, 
    add=TRUE, col="#33FF00FF", lwd = 3)
roc(bc.dt$diagnosis[-train.idx], as.vector(bc.ridge.pred.obs), plot=TRUE, legacy.axes=TRUE, 
    add=TRUE, col="#0066FFFF", lwd = 3)
roc(bc.dt$diagnosis[-train.idx], as.vector(bc.lasso.pred.obs), plot=TRUE, legacy.axes=TRUE, 
    add=TRUE, col="#FF0099FF", lwd = 3)
legend("bottomright", c("backwards elimination: AUC = 0.991", "forward selection: AUC = 0.991", "ridge regression: AUC = 0.988", "lasso regression: AUC = 0.992"), 
       fill = c("#FF9900FF","#33FF00FF","#0066FFFF", "#FF0099FF"))

# Comparing the training AUCs obtained at points (b) and (g) with the test AUCs 
# and commenting on the overfitting of each model

# predict the outcome of lasso and ridge for the observations in the training data set
bc.ridge.pred <- predict(bc.cv.ridge, newx = x.bc.dt[train.idx,], type="response", s=bc.cv.ridge$lambda.1se)
bc.lasso.pred <- predict(bc.cv.lasso, newx = x.bc.dt[train.idx,], type="response", s=bc.cv.lasso$lambda.1se)

# compute training AUCs
roc(bc.dt$diagnosis[train.idx], bc.sel.back.pred)
roc(bc.dt$diagnosis[train.idx], bc.sel.forw.pred)
roc(bc.dt$diagnosis[train.idx], as.vector(bc.ridge.pred))
roc(bc.dt$diagnosis[train.idx], as.vector(bc.lasso.pred))

# plot training AUCc
par(mfrow=c(1,1))
roc(bc.dt$diagnosis[train.idx], bc.sel.back$fitted.values, plot=TRUE, legacy.axes=TRUE, quiet = TRUE,
    main = "Training AUCs",
    col="#FF9900FF", lwd = 3)
roc(bc.dt$diagnosis[train.idx], bc.sel.forw$fitted.values, plot=TRUE, legacy.axes=TRUE, 
    add=TRUE, col="#33FF00FF", lwd = 3)
roc(bc.dt$diagnosis[train.idx], as.vector(bc.ridge.pred), plot=TRUE, legacy.axes=TRUE, 
    add=TRUE, col="#0066FFFF", lwd = 3)
roc(bc.dt$diagnosis[train.idx], as.vector(bc.lasso.pred), plot=TRUE, legacy.axes=TRUE, 
    add=TRUE, col="#FF0099FF", lwd = 3)
legend("bottomright", 
       c("backwards elimination: AUC = 0.992", "forward selection: AUC = 0.989", "ridge regression: AUC = 0.972", "lasso regression: AUC = 0.978"), 
       fill = c("#FF9900FF","#33FF00FF","#0066FFFF", "#FF0099FF"))