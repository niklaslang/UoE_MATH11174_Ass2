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
    
    # regression model summary
    log.regr.summary <- coef(summary(log.regr))
    
    # extract and calculate values from summary
    SNP <- colnames(data)[snp]
    beta <- round(log.regr.summary[2,1], 3)
    odds.ratio <- round(exp(log.regr.summary[2,1]), 3)
    std.error <- round(log.regr.summary[2,2], 3)
    p.value <- format(log.regr.summary[2,4], scientific = FALSE)
    
    # SNP summary statistics
    snp.summary <- c(SNP, beta, odds.ratio, std.error, p.value)
    
    # add SNP summary to output table
    output <- rbind(output, snp.summary)
    
  }

  output <- data.table(output)
  
  # add colnames to output table
  colnames(output) <- c("SNP", "beta", "odds.ratio", "std.error", "p.value")
  
  if(order == TRUE){
    
    # sort output by increasing p-value
    setorder(output, p.value)
  }
  
  return(output)

}

### (c) ###

univ.glm.test(gdm.dt[,!c("ID","sex","pheno")], gdm.dt$pheno, order = TRUE)
