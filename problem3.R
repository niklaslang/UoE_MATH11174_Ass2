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

snp.allele.dt <- data.table( "rsID" = str_extract(colnames(gdm.dt)[-c(1,2,3)], rsID.regex),
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
  output <- cbind(snp.allele.dt$rsID, output)
  
  # add colnames to output table
  colnames(output) <- c("SNP","beta", "std.error", "p.value")
  
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

gdm.snp.dt <- univ.glm.test(gdm.dt[,!c("ID","sex","pheno")], gdm.dt$pheno, order = TRUE)

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


