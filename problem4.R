library(data.table)

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
plot(pca.genes$x[, 1:2], main="Projection of variables on the first 2 PCs", col = "firebrick", cex=0.8)

# percentage of variance explained by the first two components
perc.expl <- pca.genes$sdev^2 / sum(pca.genes$sdev^2) 
sum(perc.expl[1:2])

# Starting from the PCA plot just produced, devising a simple rule PC2 < -10 
# that identifies the four gene expressions that are most different from the rest and report their names
pca.genes$x[which(pca.genes$x[,"PC2"] < -10), .SD]
  