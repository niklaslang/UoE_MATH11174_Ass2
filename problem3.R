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