### (a) ###

# Reading file `GDM.raw.txt` into a data table named `gdm.dt`, 
# and storing rsID and coded allele in a separate data table (call it snp.allele.dt)

library(data.table)
gdm.dt <- fread("GDM.raw.txt")

if(!"stringr" %in% rownames(installed.packages())) {
  install.packages("stringr") # functions for making regex stuff super intuitive
}

library(stringr)

rsID.regex <- "rs\\d+"
reference.allele.regex <- "[ACGT]$"

str_extract(colnames(gdm.dt)[-c(1,2,3)], rsID.regex)
str_extract(colnames(gdm.dt)[-c(1,2,3)], reference.allele.regex)

snp.allele.dt <- data.table( "snp.name" = colnames(gdm.dt)[-c(1,2,3)],
                             "rsID" = str_extract(colnames(gdm.dt)[-c(1,2,3)], rsID.regex),
                             "reference.allele" = str_extract(colnames(gdm.dt)[-c(1,2,3)], reference.allele.regex))

# Imputing missing values in `gdm.dt`` according to SNP-wise average allele count

table(is.na(gdm.dt))

for (colnm in colnames(gdm.dt[,-1])) {
  gdm.dt[[colnm]][is.na(gdm.dt[[colnm]])] <- mean(gdm.dt[[colnm]], na.rm = TRUE)
}

table(is.na(gdm.dt))

### (b) ###

# Writing a function

# univ.glm.test <- function(x, y, order=FALSE)
# where x is a data table of SNPs, y is a binary outcome vector, 
# and order is either TRUE or FALSE: 

# the function should fit a logistic regression model for each SNP in x, 
# and return a data table containing SNP names, 
# regression coefficients, 
# odds ratios, 
# standard errors 
# and p-values. 
# If order is set to TRUE, the output data table should be ordered by increasing p-value

univ.glm.test <- function( data, outcome, order=FALSE){
  
  output <- NULL
  
  # loop over all SNPs
  
  for (snp in 1:length(data)){
    
    # assertion check
    stopifnot(length(outcome) == length(data[[snp]]))
    
    # fit logistic regression model
    log.regr <- glm(outcome ~ data[[snp]], family = "binomial")
    
    # regression model summary with beta, std.error and p.value
    log.regr.summary <- data.table(signif(coef(summary(log.regr)),3))[-1,-3] # exclude intercept and t-value
    
    # add SNP summary to output table
    output <- rbind(output, log.regr.summary)
    
  }

  # add column of SNP IDs
  output <- cbind(snp.allele.dt$snp.name, output)
  
  # add colnames to output table
  colnames(output) <- c("snp","beta", "std.error", "p.value")
  
  # compute odds ratio
  output[, odds.ratio := signif(exp(beta),3)]
  
  if(order == TRUE){
    
    # sort output by increasing p-value
    setorder(output, p.value)
  }
  
  return(output)

}

### (c) ###

# Using function `univ.glm.test()`,
# running an association study for all the SNPs in `gdm.dt` against having gestational diabetes 
# (column “pheno”). 

gdm.snp.dt <- univ.glm.test(gdm.dt[,!c("ID","sex","pheno")], gdm.dt$pheno)

# For the SNP that is most strongly associated to increased risk of gestational diabetes 
# and the one with most significant protective effect, 
# reporting the summary statistics from the GWAS 
# as well as the 95% and 99% confidence intervals on the odds ratio

# SNP most strongly associated with gestational diabetes

gdm.snp.dt[p.value == min(p.value)]

beta1 <- gdm.snp.dt[beta == max(beta), beta]
se1 <- gdm.snp.dt[beta == max(beta), std.error]

# 95% CI
round(exp(beta1 + 1.96 * se1 * c(-1, 1)), 3)

# 99% CI
round(exp(beta1 + 2.58 * se1 * c(-1, 1)), 3)

# SNP with most significant protective effect

gdm.snp.dt[odds.ratio == min(odds.ratio)]

beta2 <- gdm.snp.dt[odds.ratio == min(odds.ratio), beta]
se2 <- gdm.snp.dt[odds.ratio == min(odds.ratio), std.error]

# 95% CI
round(exp(beta2 + 1.96 * se2 * c(-1, 1)), 3)

# 99% CI
round(exp(beta2 + 2.58 * se2 * c(-1, 1)), 3)

### (d) ###

# merging the GWAS results with the table of gene names provided in file `GDM.annot.txt`

gdm.annot.dt <- fread("GDM.annot.txt")

gdm.gwas.dt <- merge(snp.allele.dt, gdm.annot.dt,
                     by.x = "rsID", by.y = "snp")

gdm.gwas.dt <- merge(gdm.gwas.dt, gdm.snp.dt,
                     by.x = "snp.name", by.y = "snp")

gdm.gwas.dt[, pos := as.numeric(pos)]

# reporting SNP name, effect allele, chromosome number and corresponding gene name
# for all SNPs that have p-value < 10−4 

hit.snp.dt <- gdm.gwas.dt[p.value < 1e-4]

hit.snp.dt[,c("snp.name","reference.allele","chrom","gene")]

# for all hit SNPs reporting all gene names that are within a 1Mb window from the SNP position on the same chromosome

# hit no.1
gdm.gwas.dt[chrom == hit.snp.dt$chrom[1]][pos >= hit.snp.dt$pos[1] - 1000000 & pos <= hit.snp.dt$pos[1] + 1000000]$gene

# hit.no.2
gdm.gwas.dt[chrom == hit.snp.dt$chrom[2]][pos >= hit.snp.dt$pos[2] - 1000000 & pos <= hit.snp.dt$pos[2] + 1000000]$gene


### (e) ###

# Building a weighted genetic risk score that includes all SNPs with p-value < 10−4,
# a second score with all SNPs with p-value < 10−3, 
# and a third score that only includes SNPs on the FTO gene 

# ensure that the ordering of SNPs is respected
gdm.gwas.dt <- gdm.gwas.dt[match(colnames(gdm.dt)[-c(1,2,3)], gdm.gwas.dt$snp.name),]

# assertion check that the ordering of SNPs is respected
stopifnot(colnames(gdm.dt)[-c(1,2,3)] == gdm.gwas.dt$snp.name)

# score 1: p.value < 10^-4
gdm1.snp <- gdm.gwas.dt[p.value < 1e-4]
gdm1.grs <- gdm.dt[, .SD, .SDcols = gdm.gwas.dt[p.value < 1e-4]$snp.name]
gdm1.weighted.grs <- as.matrix(gdm1.grs) %*% gdm1.snp$beta

# score 2: p.value < 10^-3
gdm2.snp <- gdm.gwas.dt[p.value < 1e-3]
gdm2.grs <- gdm.dt[, .SD, .SDcols = gdm.gwas.dt[p.value < 1e-3]$snp.name]
gdm2.weighted.grs <- as.matrix(gdm2.grs) %*% gdm2.snp$beta

# score 3: SNP on the FTO gene
gdm3.snp <- gdm.gwas.dt[gene == "FTO"]
gdm3.grs <- gdm.dt[, .SD, .SDcols = gdm.gwas.dt[gene == "FTO"]$snp.name]
gdm3.weighted.grs <- as.matrix(gdm3.grs) %*% gdm3.snp$beta

# adding the three scores as columns to the `gdm.dt` data table
gdm.dt$p4.score <- gdm1.weighted.grs
gdm.dt$p3.score <- gdm2.weighted.grs
gdm.dt$FTO.score <- gdm3.weighted.grs

# fitting the three scores in separate logistic regression models to test their association 
# with gestational diabetes, and for each report odds ratio, 95% confidence interval and p-value

# function to fit logistic regression models for and calculate odds ratio, 95% CI l and p-value
fit.score.model <- function(formula, data){

  # fit logistic regression model  
  log.regr <- glm(formula, data = data, family = "binomial")
  
  # compute odds ratio, 95% CI and p-value
  odds.ratio <- exp(coef(summary(log.regr))[2,1])
  CI.lower <- exp(confint(log.regr)[2,1])
  CI.upper <- exp(confint(log.regr)[2,2])
  p.value <- coef(summary(log.regr))[2,4]
  
  # brief summary table of the summary statistics
  gdm.grs.dt <- data.table(rbind(NULL, c(round(odds.ratio,3), round(CI.lower,3), round(CI.upper,3), signif(p.value,3))))
  colnames(gdm.grs.dt) <- c("odds.ratio","2.5%","97.5%","p.value")
  
  return(gdm.grs.dt)
}

# fitting the three scores in separate logistic regression models
# score 1: SNPs with p.value < 10^-4
fit.score.model(pheno ~ p4.score, gdm.dt)
# score 2: SNPs with p.value < 10^-3
fit.score.model(pheno ~ p3.score, gdm.dt)
# score 3: SNPs on the FTO gene
fit.score.model(pheno ~ FTO.score, gdm.dt)

### (f) ###

# Reading the file `GDM.test.txt` into variable `gdm.test.dt`

gdm.test.dt <- fread("GDM.test.txt", stringsAsFactors = TRUE)

# For the set of patients in `gdm.test.dt`,
# computing the three genetic risk scores as defined at point (e) using the same set of SNPs and corresponding weights

# ensure that the ordering of SNPs is respected
gdm.gwas.dt <- gdm.gwas.dt[match(colnames(gdm.test.dt)[-c(1,2,3)], gdm.gwas.dt$rsID),]

# assertion check that the ordering of SNPs is respected
stopifnot(colnames(gdm.test.dt)[-c(1,2,3)] == gdm.gwas.dt$rsID)

# score 1: p.value < 10^-4
gdm.test1.grs <- gdm.test.dt[, .SD, .SDcols = gdm.gwas.dt[p.value < 1e-4]$rsID]
gdm.test1.weighted.grs <- as.matrix(gdm.test1.grs) %*% gdm1.snp$beta

# score 2: p.value < 10^-3
gdm.test2.grs <- gdm.test.dt[, .SD, .SDcols = gdm.gwas.dt[p.value < 1e-3]$rsID]
gdm.test2.weighted.grs <- as.matrix(gdm.test2.grs) %*% gdm2.snp$beta

# score 3: SNP on the FTO gene
gdm.test3.grs <- gdm.test.dt[, .SD, .SDcols = gdm.gwas.dt[gene == "FTO"]$rsID]
gdm.test3.weighted.grs <- as.matrix(gdm.test3.grs) %*% gdm3.snp$beta

# Adding the three scores as columns to `gdm.test.dt` (hint: use the same column names as before)

gdm.test.dt$p4.score <- gdm.test1.weighted.grs
gdm.test.dt$p3.score <- gdm.test2.weighted.grs
gdm.test.dt$FTO.score <- gdm.test3.weighted.grs

### (g) ###



