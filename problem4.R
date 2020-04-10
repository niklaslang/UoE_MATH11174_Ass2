library(data.table)

### (a) ###

# Reading file `nki.csv` into a data table named `nki.dt`:

nki.dt <- data.table(read.csv("nki.csv", stringsAsFactors = TRUE))

# Computing the matrix of correlations between the gene expression variables, 
# and displaying it so that a block structure is highlighted:

# library for creating corrplots
library(corrplot)

# computing the correlation between the numerical colums `nki.dt`
nki.corr <- cor(nki.dt[, .SD, .SDcols = sapply(nki.dt, is.numeric)], use="pairwise.complete")

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
        distinct.unique.pair <- c(rownames(nki.corr)[i], colnames(nki.corr)[j],nki.corr[i,j])
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

