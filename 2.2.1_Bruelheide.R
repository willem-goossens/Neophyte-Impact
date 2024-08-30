

setwd('C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/EIVE Data/ESy-master')
source('code/prep.R') #### Loading packages

### define expert file:
expertfile <- "EUNIS-ESy-2020-06-08.txt" # latest version of the EUNIS system, used in Chytrý et al. 2020, AVS

### Start  ######################
### Read and parse the expert file
#################################################################### #
### Step 1: Parse membership formulas into membership expressions ####
### Step 2: Add right-hand sides of membership expressions        ####
###         where there are no right-hand side conditions
#################################################################### #
source('code/step1and2_load-and-parse-the-expert-file.R')

######################################################## #
### Input 2                                           ####
### Read vegetation data file into a long table with class data.table and with columns named
### "RELEVE_NR", "TaxonName", "Cover_Perc"
### If you want to use your own data, please make sure to have the three columns.




###################################### #
### Step 4: Aggregate taxon levels  ####
###################################### #
# also changed file to return all non found species as na
source('code/step4_aggregate-taxon-levels.R')

setwd("C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/Neophyte-Impact")
stopCluster(cl)
