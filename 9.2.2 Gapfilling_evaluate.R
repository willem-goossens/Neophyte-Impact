
# This file must be run in R 3.4.4
library(dplyr)
library(tibble)
library(readr)
library(BHPMF)
  
######Basic input data
trait <- read_csv("../Extra data/Traits/Eva_diaz.csv")
trait <- trait[!is.na(trait$Genus),]
#trait <- trait[sample(nrow(trait), round(nrow(trait)/50), replace = F),]

# change wrong family
trait$Family[trait$Genus=="Antirrhinum"] <- gsub("Scrophulariaceae","Plantaginaceae", trait$Family[trait$Genus=="Antirrhinum"])
colnames(trait)

# create hierarchy
hier <- trait[, c(1, 17,16)]
colnames(hier)[1] <- "SpeciesName"
hier$IDX <- c(1: nrow(hier))
hier <- hier[, c(4,1,2,3)]
phylo <- read_csv("../Extra data/Species names/phylo.csv")
duplicates <- unique(hier[, 3:4])
duplicates <- duplicates[duplicated(duplicates$Genus)|duplicated(duplicates$Genus,fromLast=T),]

change <- data.frame(Genus = c("Holboellia","Calophyllum","Labisia","Platysace","Ripogonum","Sarcotheca","Sphenostemon","Stauntonia","Tectaria","Celtis","Digitalis","Hypericum","Melampyrum"),
                     Family= c("Lardizabalaceae","Calophyllaceae","Primulaceae","Apiaceae","Ripogonaceae","Oxalidaceae","Paracryphiaceae","Lardizabalaceae","Polypodiaceae","Ulmaceae","Plantaginaceae","Hypericaceae","Orobanchaceae"))

hier$Family[hier$Genus %in% change$Genus] <- change$Family[match(hier$Genus[hier$Genus %in% change$Genus], change$Genus)]
colnames(trait)
hier <- as.matrix(hier)
trait <- as.matrix(trait[, c(7:12)])

back_trans_pars <- list()
rm_col <- c()
for(i in 1:ncol(trait)){
  x <- trait[,i] # goes through the columns
  min_x <- min(x,na.rm = T) # takes the min of each column
  if(min_x < 0.00000000001){
    x <- x - min_x + 1 # make this optional if min x is neg
  }
  logx <- log10(x)
  mlogx <- mean(logx, na.rm = T)
  slogx <- sd(logx, na.rm = T)
  x <- (logx - mlogx)/slogx # Z transformation
  back_trans_pars[[i]] <- list(min_x = min_x,
                               mlogx = mlogx,
                               slogx = slogx)
  trait[,i] <- x
}

#saveRDS(back_trans_pars, "../Extra data/Traits/back_trans_pars.rds")


###########10-fold imputation to calculate RMSE
####BHPMF output RMSE automatically
#BHPMF for 4397 species
samplenum <- 1000
burnnum <- 200
gapnum <- 2
foldnum <- 10

GapFilling(as.matrix(trait), hierarchy.info = as.matrix(hier),
           prediction.level = 4, used.num.hierarchy.levels = 3,
           num.samples=samplenum, burn=burnnum, gaps=gapnum,
           num.latent.feats=20, tuning=TRUE, num.folds.tuning=foldnum,rmse.plot.test.data=TRUE,
           mean.gap.filled.output.path = "C:/Users/u0166342/Documents/Doctoraat/Impact/Extra data/Traits/Results/mean_gap_filled.txt",
           std.gap.filled.output.path= "C:/Users/u0166342/Documents/Doctoraat/Impact/Extra data/Traits/Results/std_gap_filled.txt")

rmse_out <- CalculateCvRmse(as.matrix(trait), as.matrix(hier), 
                            num.samples=samplenum, burn=burnnum, 
                            gaps=gapnum, num.folds=10, num.latent.feats=15,
                            tmp.dir="C:/Users/u0166342/Documents/Doctoraat/Impact/Extra data/Traits/Results")

rmse_out[["avg.rmse"]]
0.4177873
rmse_out[["std.rmse"]]
0.02759093


# Assign
mean <- read.delim("../Extra data/Traits/Results/mean_gap_filled.txt")
std <- read.csv("../Extra data/Traits/Results/std_gap_filled.txt")
back_trans_pars <- readRDS( "../Extra data/Traits/back_trans_pars.rds")

# back transform
for(i in 1:6){
  mean[,i] <- 10^(mean[,i]*back_trans_pars[[i]]$slogx+back_trans_pars[[i]]$mlogx)
}

boxplot(mean[,1])
boxplot(mean[,2])
boxplot(mean[,3])
boxplot(mean[,4])
boxplot(mean[,5])
boxplot(mean[,6])

hier <- as.data.frame(hier)
mean <- cbind(hier[, -1], mean)

write_csv(mean, "../Results/trait_normalized.csv")

for(i in 1:6){
  Q <- quantile(mean[,i], probs= c(.25, .75), na.rm=T)
  iqr <- IQR(mean[,i])
  #save_mean <- mean[!is.na(trait[, 6+i]), i]
  mean[,i] <- ifelse(mean[,i] <= (Q[1] - 1.5*iqr), NA, mean[,i])
  mean[,i] <- ifelse(mean[,i] >= (Q[2] + 1.5*iqr), NA, mean[,i])
  #mean[!is.na(trait[, 6+i]), i] <- save_mean
  print(sum(is.na(mean[,i])))
}