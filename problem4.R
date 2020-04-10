library(data.table)
library(glmnet)
library(caret)

### (a) ###

# Reading file `nki.csv` into a data table named `nki.dt`:

nki.dt <- data.table(read.csv("nki.csv", stringsAsFactors = TRUE))
nki.genes.dt <- nki.dt[, -c(1,2,3,4,5,6)]

# Computing the matrix of correlations between the gene expression variables, 
# and displaying it so that a block structure is highlighted:

# library for creating corrplots
library(corrplot)

# computing the correlation between the numerical colums `nki.dt`
nki.corr <- cor(nki.genes.dt, use="pairwise.complete")

# correlation plot
par(cex=0.5)
corrplot(nki.corr, order="hclust", diag=FALSE, tl.col="black",
         title="Correlation matrix (ordered by hierarchical clustering)",
         mar=c(0,2,4,0))
#corrplot(nki.corr, order="AOE", diag=FALSE, tl.col="black")
#corrplot(nki.corr, order="FPC", diag=FALSE, tl.col="black")

# Writing some code to identify the unique pairs of (distinct) variables 
# that have correlation coefficient greater than 0.80 in absolute value 
# and reporting their correlation coefficients:

# unique: ony FTO x BRCA, not BRCA x FTO in addition
# distinct: not FTO x FTO

corr.pairs <- NULL
for (i in 1:nrow(nki.corr)){
  for (j in 1:ncol(nki.corr)){
    if(i < j){
      if(abs(nki.corr[i,j]) > 0.8){
        distinct.unique.pair <- c(rownames(nki.corr)[i], colnames(nki.corr)[j], signif(nki.corr[i,j],3))
        corr.pairs <- rbind(corr.pairs, distinct.unique.pair)
      }
    }
  }
}

corr.pairs.dt <- data.table(corr.pairs)
colnames(corr.pairs.dt) <- c("gene1", "gene2", "correlation.coeff")
setorder(corr.pairs.dt, -correlation.coeff)
corr.pairs.dt

### (b) ###

# Running PCA (only over the columns containing gene expressions) 
# so that it is possible to identify variable clusters
# Producing a scatter plot of the projection of the predictors on the first two principal components 
# and reporting the percentage of variance explained by the first two components. 

# assertion check whether there are any missing values in nki.dt
stopifnot(sum(is.na(nki.genes.dt)) == 0)

pca.genes <- prcomp(t(nki.genes.dt), scale = TRUE)

# scatter plot: projection of the predictors on the first two principal components
par(cex=1)
plot(pca.genes$x[, 1:2], main="Projection of variables on the first 2 PCs", col = "firebrick", pch=19)

# percentage of variance explained by the first two components
perc.expl <- pca.genes$sdev^2 / sum(pca.genes$sdev^2) 
sum(perc.expl[1:2])

# Starting from the PCA plot just produced, devising a simple rule PC2 < -10 
# that identifies the four gene expressions that are most different from the rest and report their names
pca.genes$x[which(pca.genes$x[,"PC2"] < -10), .SD]

### (c) ###

# Running PCA (only over the columns containing gene expressions) again,
# this time in order to derive a patient-wise summary of all gene expressions (dimensionality reduction), 
# and only keeping the first 3 principal components. 
pca.patients <- prcomp(nki.genes.dt, scale = TRUE)
plot(pca.patients$x[, 1:2], main="Projection of patients on the first 2 PCs", 
     col = ifelse(nki.dt$Event == 1, "firebrick", "royalblue1"), pch=19, lwd=0.5)
legend("bottomleft", legend=c("No complications", "Complications"), col = c("royalblue1", "firebrick"), pch = 21)

top3.PCs <- pca.patients$x[, 1:3]

# add top 3 PCs to nki.dt
nki.dt$PC1 <- top3.PCs[,1]
nki.dt$PC2 <- top3.PCs[,2]
nki.dt$PC3 <- top3.PCs[,3]

# Testing if those principal components (independently) are associated with the outcome 
# in unadjusted logistic regression models and in models adjusted for age, estrogen receptor and grade

# testing PC1
PC1.unadjusted.model <- glm(Event ~ PC1, data = nki.dt, family="binomial") # unadjusted model
PC1.adjusted.model <- glm(Event ~ PC1 + Age + EstrogenReceptor + Grade, data = nki.dt, family="binomial") # adjusted model
coef(summary(PC1.unadjusted.model))[-1,c(1,4)]
coef(summary(PC1.adjusted.model))[-1,c(1,4)]

# testing PC2
PC2.unadjusted.model <- glm(Event ~ PC2, data = nki.dt, family="binomial")  # unadjusted model
PC2.adjusted.model <- glm(Event ~ PC2 + Age + EstrogenReceptor + Grade, data = nki.dt, family="binomial")  # adjusted model
coef(summary(PC2.unadjusted.model))[-1,c(1,4)]
coef(summary(PC2.adjusted.model))[-1,c(1,4)]

# testing PC3
PC3.unadjusted.model <- glm(Event ~ PC3, data = nki.dt, family="binomial")   # unadjusted model
PC3.adjusted.model <- glm(Event ~ PC3 + Age + EstrogenReceptor + Grade, data = nki.dt, family="binomial")  # adjusted model
coef(summary(PC3.unadjusted.model))[-1,c(1,4)]
coef(summary(PC3.adjusted.model))[-1,c(1,4)]

# Justify the difference in results between unadjusted and adjusted models

### (d) ###

# Fitting a lasso model to predict the binary outcome using all available predictors, 
# learning the optimal penalty parameter that maximises the AUC

prepare.glmnet <- function(data, formula=~ .){
  
  # create the design matrix to deal correctly with factor variables, 
  # without losing rows containing NAs
  old.opts <- options(na.action='na.pass')
  x <- model.matrix(formula, data)
  options(old.opts)
  
  # remove the intercept column, as glmnet will add one by default 
  x <- x[, -match("(Intercept)", colnames(x))]
  return(x)
}

# transforming the data
# predictors
x.nki.dt <- prepare.glmnet(nki.dt[,!c("PC1","PC2","PC3")], ~ . - Event) # exclude the outcome
# outcome
y.nki.dt <- nki.dt$Event
# set.seed
set.seed(1)

# fitting a lasso model and learning by cross-validation the penalty parameter λ that maximizes the AUC
# reporting run time
system.time(nki.cv.lasso1 <- cv.glmnet(x.nki.dt, y.nki.dt, family="binomial", type.measure="auc"))

# penalty parameters λ that maximizes the AUC
nki.cv.lasso1$lambda.min

# max AUC
nki.cv.lasso1$cvm[which(nki.cv.lasso1$lambda == nki.cv.lasso1$lambda.min)]

# Repeat the same procedure but this time penalising only to the gene expression variables, 
# leaving all other covariates unpenalised 
# reporting run time
system.time(nki.cv.lasso2 <- cv.glmnet(x.nki.dt, y.nki.dt, penalty.factor = c( rep(0, 5), rep(1,71)),
                          family="binomial", type.measure="auc"))

# penalty parameters λ that maximizes the AUC
nki.cv.lasso2$lambda.min
# max AUC
nki.cv.lasso2$cvm[which(nki.cv.lasso2$lambda == nki.cv.lasso2$lambda.min)]

# Can you explain the difference between the two models in terms of AUC? 