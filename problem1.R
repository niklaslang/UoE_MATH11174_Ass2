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

lipids.dt[, lipid.class := NA]
  
for(lipid in lipid.classes.dt$CE){
  
    # create regex that finds all lipid class identifiers
    
    regex <- paste0("(\\s|^)[", 
                    toupper(substr(lipid,1,1)),
                    tolower(substr(lipid,1,1)),
                    "](",
                    toupper(substr(lipid,2,nchar(lipid))),
                    "|",
                    tolower(substr(lipid,2,nchar(lipid))),
                    ")")
    
    # add new col lipid.class to lipid.dt
    lipids.dt[, lipid.class := ifelse(grepl(regex, lipids.dt$lipid.species), lipid, lipid.class)]
}

### (b) ###

# Computing the Wald test statistic for each lipid species

lipids.dt[, "wald.test" := round((oddsratio-1)/se,3)]

# Assuming that there were 288 patients in the dataset, deriving a p-value for each lipid species 
# using the t distribution and append them to the results.dt data table

# t.distn with n=288???

lipids.dt[, "p.value(t-distn)" := signif(1-(pt(abs(wald.test), 277)-pt(-abs(wald.test), 277)),3)]

# Repeating the computation of the p-values using the normal distribution (hint: you’ll need to use the pnorm() function)

lipids.dt[, "p.value(n-distn)" := signif(1-(pnorm(abs(wald.test))-pnorm(-abs(wald.test))),3)]

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
