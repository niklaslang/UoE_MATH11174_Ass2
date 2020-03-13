### (a) ###

# load library to work with data tables
library(data.table)

# import file and look up files
lipids.dt <- fread("lipids.txt", stringsAsFactors = TRUE)
lipid.classes.dt <- fread("lipid-classes.txt", stringsAsFactors = TRUE)

# Using the lipid-class lookup file to determine the possible classes,
# then filter the first column in order to determine which class the species belongs
# then merge with the class table

# filter lipids.dt$lipid.species for lipid class

classify.lipids <- function(lipid.table, lipid.species.col, lipid.classes){
  
  # add  class colum to lipid.table
  lipid.table[, lipid.class := NA]
  
  # loop over over all lipid classes
  
  for(lipid in lipid.classes){
    
    # create regex that finds all lipid class identifiers
    regex <- paste0("(\\s|^)[", 
                    toupper(substr(lipid,1,1)),
                    tolower(substr(lipid,1,1)),
                    "](",
                    toupper(substr(lipid,2,nchar(lipid))),
                    "|",
                    tolower(substr(lipid,2,nchar(lipid))),
                    ")(\\s|$)")
  
    # add new col lipid.class to lipid.dt
    
    lipid.species <- paste0(lipid.table, "$" , lipid.species.col)
    lipid.table[, lipid.class := ifelse(grepl(regex, lipid.species), lipid, lipid.class)]
  }
}

# run function
classify.lipids(lipids.dt, "lipid.species" , lipid.classes.dt$CE)

# show results
head(lipids.dt, n = 15)

# update lipid.classes

lipid.classes.dt <- rbind(lipid.classes.dt, list("CE", "Cholesterylesters"))

# rerun fucnction
classify.lipids(lipids.dt, "lipid.species" , lipid.classes.dt$CE)

# merge data tables
results.dt <- merge(lipids.dt, lipid.classes.dt, 
                    by.x = "lipid.class", 
                    by.y = "CE" )

# Count the number of lipids that fall in each class
table(results.dt$`Cholesterol esters`)
  
# alternative: just one for loop that has to be rerun...

#lipids.dt[, lipid.class := NA]
#  
#for(lipid in lipid.classes.dt$CE){
#  
#    # create regex that finds all lipid class identifiers
#    
#    regex <- paste0("(\\s|^)[", 
#                    toupper(substr(lipid,1,1)),
#                    tolower(substr(lipid,1,1)),
#                    "](",
#                    toupper(substr(lipid,2,nchar(lipid))),
#                    "|",
#                    tolower(substr(lipid,2,nchar(lipid))),
#                    ")")
#    
#    # add new col lipid.class to lipid.dt
#    lipids.dt[, lipid.class := ifelse(grepl(regex, lipids.dt$lipid.species), lipid, lipid.class)]
#}

### (b) ###

# Computing the Wald test statistic for each lipid species

results.dt[, "wald.test" := round((oddsratio-1)/se,3)]

# Assuming that there were 288 patients in the dataset, deriving a p-value for each lipid species 
# using the t distribution and append them to the results.dt data table

# t.distn with n=288???

results.dt[, "p.value(t-distn)" := signif(1-(pt(abs(wald.test), 277)-pt(-abs(wald.test), 277)),3)]

# Repeating the computation of the p-values using the normal distribution (hint: you’ll need to use the pnorm() function)

results.dt[, "p.value(n-distn)" := signif(1-(pnorm(abs(wald.test))-pnorm(-abs(wald.test))),3)]

# Providing some evidence to justify if the normal approximation is acceptable in this instance

hist(lipids.dt$wald.test, 
     freq=F,
     breaks = 16,
     main = "Distribution of Wald test statistic",
     xlab = "Wald test statistic",
     ylab = "Frequency")

lines(density(lipids.dt$wald.test), col="red")

lines(seq(-8, 8, by=.05), 
      dnorm(seq(-8, 8, by=.05), mean(lipids.dt$wald.test), sd(lipids.dt$wald.test)),
      col="blue")

legend("topright", c("observed density", "normal density"), col = c("red","blue"), lty = c(1,1))

### (c) ###

# Implementing the following function: holm.bonferroni <- function(results.dt, alpha)
# where results.dt is a data table containing the hypotheses tested and corresponding p-values,
# and alpha is the desired significance level
# the function should return the subset of results.dt that are significant 
# according to the Holm-Bonferroni method ordered by increasing p-value

lipids.results.dt <- lipids.dt[,c('lipid.species','p.value(t-distn)')]

holm.bonferroni <- function(results.dt, alpha=0.05){
  
  # order p-values
  setorder(results.dt, cols = 'p.value(t-distn)')
  
  # variables
  
  m <- length(results.dt$`p.value(t-distn)`)
  k <- 1
  
  # loop over all p-values
  for (p.k in results.dt$`p.value(t-distn)`){
    
    # reject yes/no?
    if(p.k < alpha/(m+1-k)){
      #print(k)
      #print(p.k)
      k <- k + 1
    }else{
      break
    }
  }

  # subset output data table
  fwer.results.dt <- results.dt[1:k-1,c('lipid.species','p.value(t-distn)')]
  
  # order output data table
  setorder(fwer.results.dt, cols = 'p.value(t-distn)')
  
  # return result
  return(fwer.results.dt)
}

holm.bonferroni(results.dt)

### (d) ###

# Implementing the following function: benjamini.hochberg <- function(results.dt, q)

# that implements the Benjamini-Hochberg procedure and, 
# given data table results.dt of hypotheses and p-values and a false discovery rate q, 
# returns the subset of results.dt that are significant ordered by increasing p-value.

benjamini.hochberg <- function(results.dt, q=0.05){
  
  # sort results by p-value
  setorder(results.dt, cols = 'p.value(t-distn)')
  
  # variables
  
  m <- length(results.dt$`p.value(t-distn)`)
  k <- 1
  
  # loop over all p-values
  for (p.k in results.dt$`p.value(t-distn)`){
    
    # reject yes/no?
    if(p.k < (k/m*q)){
      #print(k)
      #print(p.k)
      k <- k + 1
    }else{
      break
    }
  }
  
  # subset output data table
  fdr.results.dt <- results.dt[1:k-1,c('lipid.species','p.value(t-distn)')]
  
  # order output data table
  setorder(fdr.results.dt, cols = 'p.value(t-distn)')
  
  # return result
  return(fdr.results.dt)
}

benjamini.hochberg(results.dt, 0.05)

### (e) ###

# Produce a volcano plot using a different colour for the lipid species 
# that are significant after controlling for the family-wise error rate 
# given a nominal significance threshold α = 0.05 
# (it is not required to display names alongside the highlighted points). 
# Use a different symbol to highlight the points that are considered significant 
# according to the Benjamini-Hochberg procedure at a false discovery rate of 1%

fwer.dt <- holm.bonferroni(results.dt, 0.05)
fdr.dt <- benjamini.hochberg(results.dt, 0.01)

plot(log(results.dt$oddsratio), -log10(results.dt$`p.value(t-distn)`), 
     main="Volcano plot",
     xlab="log Odds Ratio", 
     ylab="-log10(p-value)", 
     col=ifelse(results.dt$`p.value(t-distn)` %in% fwer.dt$`p.value(t-distn)` , "#FFB900FF","#008EFFFF"),
     pch=ifelse(results.dt$`p.value(t-distn)` %in% fdr.dt$`p.value(t-distn)`, 17, 19),
     cex=ifelse((results.dt$`p.value(t-distn)` %in% fdr.dt$`p.value(t-distn)` | results.dt$`p.value(t-distn)` %in% fwer.dt$`p.value(t-distn)`), 1.2, 0.5),
     xlim=c(-1,1),
     ylim=c(0,14))

legend("topleft", c("FWER < 5%", "FWER ≥ 5%", "FDR < 1%", "FDR ≥ 1%"), 
       col = c("#FFB900FF","#008EFFFF", "grey", "grey"),
       lwd = c(3,3,1,1),
       lty = c(1,1,NA,NA),
       pch = c(NA,NA,2,1))

