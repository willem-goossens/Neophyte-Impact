
library(devtools)
library(readr)
library(dplyr)
library(BHPMF)

set.seed(123)

######Basic input data
trait <- read_csv("../TRY/Species_traits_EVA.csv", show_col_types = FALSE)
#trait <- trait[sample(nrow(trait), round(nrow(trait)/50), replace = F),]
trait$Family[trait$Genus=="Antirrhinum"] <- gsub("Scrophulariaceae","Plantaginaceae", trait$Family[trait$Genus=="Antirrhinum"])
colnames(trait)
hier <- trait[, c(1, 3,4)]
colnames(hier)[1] <- "SpeciesName"
hier$IDX <- c(1: nrow(hier))
hier <- hier |> relocate(IDX, .before= SpeciesName)
phylo <- read_csv("phylo.csv",show_col_types = FALSE)
duplicates <- unique(hier[, 3:4])
duplicates <- duplicates[duplicated(duplicates$Genus)|duplicated(duplicates$Genus,fromLast=T),]

test <- hier[hier$Genus %in% duplicates$Genus,]

test <- test |> group_by(Genus, Family) |> summarise(n=n())
test <- test[order(test$n, decreasing=T),]
test <- test[!duplicated(test$Genus),]

hier$Family[hier$Genus %in% test$Genus] <- test$Family[match(hier$Genus[hier$Genus %in% test$Genus], test$Genus)]

change <- data.frame(Genus = c("Holboellia","Calophyllum","Labisia","Platysace","Ripogonum","Sarcotheca","Sphenostemon","Stauntonia","Tectaria"),
                     Family= c("Lardizabalaceae","Calophyllaceae","Primulaceae","Apiaceae","Ripogonaceae","Oxalidaceae","Paracryphiaceae","Lardizabalaceae","Polypodiaceae"))

hier$Family[hier$Genus %in% change$Genus] <- change$Family[match(hier$Genus[hier$Genus %in% change$Genus], change$Genus)]

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
#write.csv(back_trans_pars, paste0("back_trans_pars.csv"))

###########10-fold imputation to calculate RMSE
####BHPMF output RMSE automatically
#BHPMF for 4397 species
samplenum <- 1000
burnnum <- 200
gapnum <- 2
foldnum <- 10
GapFilling(as.matrix(trait), hierarchy.info = hier,
           prediction.level = 4, used.num.hierarchy.levels = 3,
           num.samples=samplenum, burn=burnnum, gaps=gapnum,
           num.latent.feats=20, tuning=TRUE, num.folds.tuning=foldnum,rmse.plot.test.data=TRUE,
           mean.gap.filled.output.path = "C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/Intermediate Data/mean_gap_filled.txt",
           std.gap.filled.output.path= "C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/Intermediate Data/std_gap_filled.txt")

rmse_out <- CalculateCvRmse(as.matrix(trait), hier, 
                            num.samples=samplenum, burn=burnnum, 
                            gaps=gapnum, num.folds=10, num.latent.feats=15,
                            tmp.dir="C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/Intermediate Data")


# Assign
mean <- read.delim("C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/Intermediate Data/mean_gap_filled.txt")
#back_trans_pars <- read.csv("back_trans_pars.csv")

for(i in 1:6){
  mean[,i] <- 10^(mean[,i]*back_trans_pars[[i]]$slogx+back_trans_pars[[i]]$mlogx)
}



boxplot(mean[,1])
boxplot(mean[,2])
boxplot(mean[,3])
boxplot(mean[,4])
boxplot(mean[,5])
boxplot(mean[,6])



for(i in 1:6){
  Q <- quantile(mean[,i], probs= c(.25, .75), na.rm=T)
  iqr <- IQR(mean[,i])
  #save_mean <- mean[!is.na(trait[, 6+i]), i]
  mean[,i] <- ifelse(mean[,i] <= (Q[1] - 1.5*iqr), NA, mean[,i])
  mean[,i] <- ifelse(mean[,i] >= (Q[2] + 1.5*iqr), NA, mean[,i])
  #mean[!is.na(trait[, 6+i]), i] <- save_mean
  print(sum(is.na(mean[,i])))
}

mean <- cbind(hier[, -1], mean)

write_csv(mean, "trait_normalized.csv")
