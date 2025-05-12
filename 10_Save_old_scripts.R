
## 3.4 Genus traits
```{r}
# remove also aggr. from species
eva$Genus <- gsub("([A-Za-z]+).*", "\\1", eva$name)
eva$new_genus <- gsub("([A-Za-z]+).*", "\\1", eva$new_species)
eva$new_genus <- as.character(eva$new_genus)

eva_names <- unique(eva$species)

eva$official <- ifelse(!is.na(eva$new_species), eva$new_species, eva$species)
eva$official_genus <- ifelse(!is.na(eva$new_species), eva$new_genus, eva$genus)

# we get Diaz genus level data
Diaz_genus <- Diaz_raw |> group_by(Genus) |> summarise(LA= mean(`Leaf area (mm2)`, na.rm=T),
                                                       N= mean(`Nmass (mg/g)`, na.rm=T),
                                                       LMA= mean(`LMA (g/m2)`, na.rm=T), 
                                                       H= mean(`Plant height (m)`, na.rm=T),
                                                       SM= mean(`Diaspore mass (mg)`, na.rm=T),
                                                       SSD= mean(`SSD combined (mg/mm3)`, na.rm=T))


# we make a copy of eva to use as a working version
copy<- eva

# left join with Diaz
# copy <- left_join(copy, Diaz[,c(2,4,6:14)], by= c("species"="Species name standardized against TPL"))

# left join with diaz raw by official name
diaz_eva<-  full_join(eva, Diaz_raw, by= c("official"="Species name standardized against TPL"))

# change genus for species that are not present in eva
diaz_eva$official_genus[is.na(diaz_eva$official_genus)]<-
  diaz_eva$Genus[is.na(diaz_eva$official_genus)]
# give family to theses species
diaz_eva$family[is.na(diaz_eva$family)]<-
  diaz_eva$Family[is.na(diaz_eva$family)]


old_method=F
if(old_method){
  diaz_eva<- full_join(eva, Diaz_raw, by= c("new_species"="Species name standardized against TPL"))
  
  # add traits from old names not in accepted names
  diaz_eva[is.na(diaz_eva$`Number of traits with values`), c(27:39)] <-
    Diaz_raw[match(diaz_eva$species[is.na(diaz_eva$`Number of traits with values`)], 
                   Diaz$`Species name standardized against TPL`[is.na(diaz_eva$`Number of traits with values`)]), c(1,3:14)]
  
  # get all diaz eva names
  diaz_eva_names <- unique(diaz_eva$species[!is.na(diaz_eva$`Number of traits with values`)])
}

copy <- diaz_eva
copy <- copy[!is.na(copy$`Number of traits with values`),]
copy <- copy[!(copy$`Number of traits with values`==1 & !is.na(copy$`SSD combined (mg/mm3)`)),]
copy <- copy[,-c(38)]
copy <- copy[!duplicated(copy$official),]

# create list of hierarchy in eva
hierarchy.info<- copy[, c("official","official_genus","family")]
hierarchy.info$id <- c(1: nrow(hierarchy.info))
hierarchy.info <- hierarchy.info %>% relocate(id)
hierarchy.info <- hierarchy.info[!duplicated(hierarchy.info$official),]
colnames(hierarchy.info)[2] <- "species"
colnames(hierarchy.info)[3] <- "genus"
hierarchy.info<- as.data.frame(hierarchy.info)

# get trait data
trait.info <- copy[, c(33:38)]
trait.info <- as.matrix(trait.info)

trait.info<- trait.info[,c(1:5)]

back_trans_pars <- list()
rm_col <- c()
for(i in 1:ncol(trait.info)){
  x <- trait.info[,i] # goes through the columns
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
  trait.info[,i] <- x
}
```


```{r}
unique_geni <- unique(hierarchy.info$genus)

x<- data.frame(genus= character(), n= numeric(), family1= character(), n1=numeric(), family2= character(), n2=numeric(), rel1= numeric())
for(i in unique_geni){
  hier<- hierarchy.info[hierarchy.info$genus==i,]
  length <- length(unique(hierarchy.info$family[hierarchy.info$genus==i]))
  if(length > 1){
    test <- hier |> group_by(family) |> summarise(n= n()) 
    test <- test[order(-test$n),]
    x <- add_row(x, genus= i, n= length, family1= test$family[1],n1= test$n[1], family2= test$family[2], n2= test$n[2], rel1= n1/(n1+n2))
  }
}

family<- x[, c(1,3)]
colnames(family)[2]<-"family"


own_family <- data.frame(genus=c("Calophyllum","Holboellia","Ripogonum","Tectaria","Labisia",
                                 "Sarcotheca","Sphenostemon","Stauntonia"),
                         family=c("Calophyllaceae","Lardizabalaceae","Rhipogonaceae","Polypodiaceae",
                                  "Primulaceae","Anacardiaceae","Paracryphiaceae","Lardizabalaceae"))

family$family[family$genus %in% own_family$genus] <- 
  own_family$family[match(family$genus[family$genus %in% own_family$genus], own_family$genus)]

hierarchy.info$family[hierarchy.info$genus %in% family$genus]<- 
  family$family[match(hierarchy.info$genus[hierarchy.info$genus %in% family$genus],family$genus)]
```




```{r}
GapFilling(trait.info, hierarchy.info,
           mean.gap.filled.output.path= "../TRY/Species_mean_traits.txt",
           std.gap.filled.output.path="../TRY/Species_std_traits.txt")

traits <- read.delim("../TRY/Species_mean_traits.txt")
colnames(traits)[1:5]<- c("LA", "N", "LMA","H","SM")
traits <- cbind(hierarchy.info, traits)
```

# check names
eva$eive_name <- eive$TaxonConcept[match(eva$name, eive$TaxonConcept)]
new_names <- eva[is.na(eva$eive_name),]




# which species are additionally present when using species --> 345 additional species
extra<- new_names[!is.na(match(new_names$`species`, eive$TaxonConcept)),]
extra <- unique(extra[,2:6])
extra_species <- extra

# Second we add species data for these species
eva$EIVEresM[eva$species %in% extra$species] <- eive$`EIVEres-M`[match(eva$`species`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresN[eva$species %in% extra$species] <- eive$`EIVEres-N`[match(eva$`species`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresR[eva$species %in% extra$species] <- eive$`EIVEres-R`[match(eva$`species`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresL[eva$species %in% extra$species] <- eive$`EIVEres-L`[match(eva$`species`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresT[eva$species %in% extra$species] <- eive$`EIVEres-T`[match(eva$`species`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]

# niche widths
eva$EIVEnwM[eva$species %in% extra$species] <- eive$`EIVEres-M.nw3`[match(eva$`species`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwN[eva$species %in% extra$species] <- eive$`EIVEres-N.nw3`[match(eva$`species`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwR[eva$species %in% extra$species] <- eive$`EIVEres-R.nw3`[match(eva$`species`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwL[eva$species %in% extra$species] <- eive$`EIVEres-L.nw3`[match(eva$`species`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwT[eva$species %in% extra$species] <- eive$`EIVEres-T.nw3`[match(eva$`species`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]

# check new names
eva$eive_name[eva$species %in% extra$species]  <- eive$TaxonConcept[match(eva$species[eva$species %in% extra$species] , eive$TaxonConcept)]
length(eive$TaxonConcept[match(eva$irena[eva$species %in% extra_species$species] , eive$TaxonConcept)])
new_names <- eva[is.na(eva$eive_name),]
# check previous annotation
extra<- new_names[!is.na(match(new_names$`species`, eive$TaxonConcept)),]




# which species are additionally present when using irena --> 132 additional species
extra<- new_names[!is.na(match(new_names$`irena`, eive$TaxonConcept)),]
extra <- unique(extra[,2:6])
extra_irena <- extra

# Third we add eva irena data for these species
eva$EIVEresM[eva$species %in% extra$species] <- eive$`EIVEres-M`[match(eva$`irena`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresN[eva$species %in% extra$species] <- eive$`EIVEres-N`[match(eva$`irena`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresR[eva$species %in% extra$species] <- eive$`EIVEres-R`[match(eva$`irena`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresL[eva$species %in% extra$species] <- eive$`EIVEres-L`[match(eva$`irena`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresT[eva$species %in% extra$species] <- eive$`EIVEres-T`[match(eva$`irena`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]

# niche widths
eva$EIVEnwM[eva$species %in% extra$species] <- eive$`EIVEres-M.nw3`[match(eva$`irena`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwN[eva$species %in% extra$species] <- eive$`EIVEres-N.nw3`[match(eva$`irena`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwR[eva$species %in% extra$species] <- eive$`EIVEres-R.nw3`[match(eva$`irena`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwL[eva$species %in% extra$species] <- eive$`EIVEres-L.nw3`[match(eva$`irena`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwT[eva$species %in% extra$species] <- eive$`EIVEres-T.nw3`[match(eva$`irena`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
# check new names
eva$eive_name[eva$species %in% extra$species]  <- eive$TaxonConcept[match(eva$`irena`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
length(eive$TaxonConcept[match(eva$`irena`[eva$species %in% extra$species],eive$TaxonConcept)])
new_names <- eva[is.na(eva$eive_name),]
# check previous annotation
extra<- new_names[!is.na(match(new_names$`irena`, eive$TaxonConcept)),]



# finally lets check the 'Matched concept' name --> 6 additional species
extra<- new_names[!is.na(match(new_names$`Matched concept`, eive$TaxonConcept)),]
extra <- unique(extra[,2:4])
extra_matched <- extra

# Third we add eva matched concept data for these species
eva$EIVEresM[eva$species %in% extra$species] <- eive$`EIVEres-M`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresN[eva$species %in% extra$species] <- eive$`EIVEres-N`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresR[eva$species %in% extra$species] <- eive$`EIVEres-R`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresL[eva$species %in% extra$species] <- eive$`EIVEres-L`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]
eva$EIVEresT[eva$species %in% extra$species] <- eive$`EIVEres-T`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                       eive$TaxonConcept)]

# niche widths
eva$EIVEnwM[eva$species %in% extra$species] <- eive$`EIVEres-M.nw3`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwN[eva$species %in% extra$species] <- eive$`EIVEres-N.nw3`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwR[eva$species %in% extra$species] <- eive$`EIVEres-R.nw3`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwL[eva$species %in% extra$species] <- eive$`EIVEres-L.nw3`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
eva$EIVEnwT[eva$species %in% extra$species] <- eive$`EIVEres-T.nw3`[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
# check new names
eva$eive_name[eva$species %in% extra$species]  <- eive$TaxonConcept[match(eva$`Matched concept`[eva$species %in% extra$species],
                                                                          eive$TaxonConcept)]
length(eive$TaxonConcept[match(eva$`Matched concept`[eva$species %in% extra$species],eive$TaxonConcept)])
new_names <- eva[is.na(eva$eive_name),]

# check previous annotation
extra<- new_names[!is.na(match(new_names$`Matched concept`, eive$TaxonConcept)),]















hist(speciesDominance$numberOfPlots[speciesDominance$numberOfPlots<50])
sum(speciesDominance$numberOfPlots>50)
sum(speciesDominance$numberOfPlots>10)

x <- data.frame(c(10,20,30,40,50))
colnames(x) <-"length"
for(i in 1:5) {
  x$rel[i] <- sum(speciesDominance$numberOfPlots >= x$length[i])/length(speciesDominance$numberOfPlots)
}

p <-ggplot(x, aes(x= length, y= rel, label= round(rel, 2)))+
  geom_line()+
  geom_point()+
  labs(y="relative number of plots selected", x="minimum number of plots")+
  geom_label(aes(length,rel))
p
writeClipboard(p)

plot(x, type= "l")




joinedData2 <- st_join(plotLocations, eco, join = st_within)

#------------------------------------------------------------------------------#
###### 5.2 unassigned ######
#------------------------------------------------------------------------------#
percentNotAssignedPlots <- sum(is.na(joinedData$ECO_NAME))/length(joinedData$PlotObservationID)*100
percentNotAssignedPlots


# Prepare parallel
parallel::detectCores()
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

# Create threshold to use
distanceThreshold <- set_units(10000,m)
# Dataset of points that were not assigned with sf_within
remainingPlots <-plotLocations[is.na(joinedData$ECO_NAME),]

# make smaller to make calculation quicker
remainingPlots<- remainingPlots[, c("PlotObservationID","Country", "geometry")]

# Compute bounding boxes for the regions, this allows to decide faster, if a plot is outside the distanceThreshold from a region
boundingBoxes <- eco

# Create box for every country
for(i in 1:length(boundingBoxes$ECO_NAME)) {
  boundingBoxes$geometry[i] <- st_as_sfc(st_bbox(boundingBoxes$geometry[i]))
}

# Begin parameters
remainingPlots$eco <- NA
remainingPlots$eco_N <- NA
remainingPlots$Distance <- -1

# Start clock
begin<-Sys.time() 

# Perform parallel loop and add results to remaining plots
remainingPlots[,4:6]<-foreach(i = 1:length(remainingPlots$PlotObservationID), 
                              .combine='rbind', .packages=c("dplyr","mgcv", "sf","units")) %dopar% {
                                
                                # set region to no region and distance to large value
                                region <- NA
                                eco_N <-NA
                                distance <- set_units(10000000,m)
                                
                                # check for every country individually if distance to the box is more than is allowed
                                for(j in 1:length(boundingBoxes$ECO_NAME)) {
                                  sf_use_s2(F)
                                  distBoundingBox <- st_distance(boundingBoxes[j,], remainingPlots[i,]) #Calculate the distance
                                  # Only if the distance to the bounding box is smaller than distanceThreshold we need to compute the real distance to the region
                                  if(distBoundingBox < distanceThreshold){
                                    dist <- st_distance(eco[j,], remainingPlots[i,])
                                    # if this distance is smaller than the previously computed distance and the threshold, save this one
                                    if(dist < distanceThreshold & dist < distance) {
                                      region <- eco$ECO_NAME[j]
                                      eco_N <- eco$ECO_ID
                                      distance <- dist
                                    }
                                  }
                                }
                                # add these to the remaining plots
                                remainingPlots$eco[i] <- c(region, eco_N, distance)
                              }
# End time and calculation time needed
end<-Sys.time()
round(end-begin, 2)  

st_distance(eco[j,], remainingPlots[91,])

# Stop cluster
parallel::stopCluster(cl = my.cluster)

joinedData$ECO_NAME[is.na(joinedData$ECO_NAME)] <- remainingPlots$eco



#### 
# Helper function to add rows if the index is within bounds
add_row_if_exists <- function(res, baseModel, idx, class, coverClass, plantName, plantStatus, tmpFullPlotData, plotsWherePlantOccurs) {
  if (idx <= nrow(baseModel$p.table)) {
    res <- add_row(res, 
                   taxa = plantName, 
                   Estimate = baseModel$p.table[idx, 1],
                   StdErr = baseModel$p.table[idx, 2],
                   zValue = baseModel$p.table[idx, 3], 
                   pr = baseModel$p.table[idx, 4],
                   Intercept = baseModel$p.table[1, 1],  
                   numberOfPlots = length(plotsWherePlantOccurs$PlotObservationID), 
                   Neophyte = plantStatus,
                   RelDiff = 
                     (mean(fitted(base)[tmpFullPlotData$coverClass == coverClass]) -
                        mean(fitted(base)[tmpFullPlotData$coverClass == "0%"])) /
                     mean(fitted(base)[tmpFullPlotData$coverClass == "0%"]),
                   class = class, n = sum(tmpFullPlotData$coverClass == coverClass),
                   size= size, rel= length(plotsWherePlantOccurs$PlotObservationID)/ size)
  }
  return(res)
}



####Old visualisation ####



```{r}

obj <- data.frame(
  cover_rel_intra = median(fullPlotData$cover_rel_intra),
  cover_rel_extra = median(fullPlotData$cover_rel_extra),
  hfp = median(fullPlotData$hfp),
  EIVEresM= median(fullPlotData$EIVEresM),
  EIVEresN= median(fullPlotData$EIVEresN),
  EIVEresL= median(fullPlotData$EIVEresL),
  EIVEresR= median(fullPlotData$EIVEresR),
  EIVEresT= median(fullPlotData$EIVEresT),
  DistSeverity.sqrt =  median(fullPlotData$DistSeverity.sqrt),
  Soil.Disturbance.sqrt= median(fullPlotData$Soil.Disturbance.sqrt),
  Grazing.Pressure.sqrt= median(fullPlotData$Grazing.Pressure.sqrt),
  Mowing.Frequency.sqrt= median(fullPlotData$Mowing.Frequency.sqrt),
  chelsaP= median(fullPlotData$chelsaP),
  Longitude= median(fullPlotData$Longitude),
  Latitude= median(fullPlotData$Latitude),
  Dataset= names(which.max(table(fullPlotData$Dataset)))
)

test <-summary(MDL)
index <- which(unique(rownames(test$coefficients[[1]]))=="cover_rel_intra")

int <-predict(MDL, obj, type="response")
slope <- test$coefficients$cond[index]

impact_intra <- slope/int

```




Visualise
```{r}
pred1=F
if(pred1){
  # predict responses and plot
  pr <- predict_response(MDL, "cover_rel_intra [all]")
  plot(pr, show_data = T, show_residuals = T, show_residuals_line = T)
  pr <- predict_response(MDL, "cover_rel_extra [all]")
  plot(pr, show_data = TRUE, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "EIVEresN [all]")
  plot(pr, show_data = TRUE, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "EIVEresM [all]")
  plot(pr, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "EIVEresR [all]")
  plot(pr, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "EIVEresT [all]")
  plot(pr, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "EIVEresL [all]")
  plot(pr, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "hfp [all]")
  plot(pr, show_data = TRUE, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "Soil.Disturbance.sqrt  [all]")
  plot(pr, show_data = TRUE, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "DistSeverity.sqrt  [all]")
  plot(pr, show_data = TRUE, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "Grazing.Pressure.sqrt  [all]")
  plot(pr, show_data = TRUE, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "Mowing.Frequency.sqrt  [all]")
  plot(pr, show_data = TRUE, show_residuals = TRUE, show_residuals_line = TRUE)
  pr <- predict_response(MDL, "chelsaP [all]")
  plot(pr, show_data = T, show_residuals = TRUE, show_residuals_line = TRUE)
}

# plot using visreg
visreg::visreg(MDL, xvar="EIVEresN", gg = TRUE, partial=FALSE, rug = F,  scale="response")
visreg::visreg(MDL, xvar="cover_rel_intra", gg = TRUE, partial=FALSE, rug = F,  scale="response")
visreg::visreg(MDL, xvar="cover_rel_extra", gg = TRUE, partial=FALSE, rug = F,  scale="response")
```


#### 2.1 Function 1 ####
```{r, meassage=F}
# The name of the first species (to test the function)
plantName = allPlants$names[339]
plantStatus= allPlants$Neophyte[339]

which(allPlants$names=="Trientalis europaea")
# Function to compute the impact for each species, so we only need the plant name as input
computePlantImpact1 <- function(plantName, plantStatus) {
  
  # Copy dataset temporarily to be able to make changes
  tmpFullPlotData <- fullPlotData
  
  reducedEva <- reducedEva_official
  # Select only explanatory variables
  tmpFullPlotData <- tmpFullPlotData %>% 
    select(numberOfVascularPlantSpecies,PlotObservationID, Area, 
           EIVEresM, EIVEresN, EIVEresR, EIVEresL, EIVEresT, 
           DistSeverity.sqrt, Soil.Disturbance.sqrt, 
           Grazing.Pressure.sqrt, Mowing.Frequency.sqrt, 
           Latitude, Longitude, Dataset, totalCover, chelsaP, hfp, PC1,PC2,PC3,PC4)
  
  # Copy the full dataset without aberrant totalCover
  tmpFullPlotData <- tmpFullPlotData %>% filter(totalCover <= 900)
  reducedEva <- reducedEva[reducedEva$PlotObservationID %in% tmpFullPlotData$PlotObservationID,]
  
  # Get the indices of the plant in the eva dataset
  indexOfPlant <- reducedEva$name == plantName & reducedEva$Neophyte== plantStatus
  
  # We select only the rows containing the species first because this is much faster than doing it for all plots.
  # We make the dataset smaller (only those plots where the species is present), group the data per plot and compute the total cover     (all species) and cover of the species we are highlighting now
  plotsWherePlantOccurs <- reducedEva[indexOfPlant,c(1,4)]
  colnames(plotsWherePlantOccurs)[2]<- "speciesCover"
  
  # Join the computed cover plots with the full data
  tmpFullPlotData <- left_join(tmpFullPlotData, plotsWherePlantOccurs, 
                               by = "PlotObservationID")
  
  # Identify all plots where the species investigated occurs
  tmpFullPlotData$plantOccurs = tmpFullPlotData$PlotObservationID %in% reducedEva$PlotObservationID[indexOfPlant]
  
  
  # reduce the dataset to obtain only those plots in which the species is present
  species_plots <- tmpFullPlotData[tmpFullPlotData$plantOccurs, c(19:22)  ]
  species_plot_number <- nrow(species_plots)
  
  # make matrix and reduce duplicated columns
  species_plots<- as.matrix(species_plots)
  duplicated<- duplicated(species_plots)
  species_plots<- species_plots[!duplicated, ]
  
  # take ConvexHull around the PCA coordinates in which the species is present
  convex_hull <- convhulln((species_plots))
  # Extract all points from the dataset that are additionally within this ConvexHull and do not contain the species
  is_within_hull <- inhulln(convex_hull, as.matrix(tmpFullPlotData[, c(19:22)]))
  subset_data <- tmpFullPlotData[is_within_hull, ]
  
  # Take size of subset data
  size<- length(subset_data$plantOccurs)
  species_plot_number <- sum(subset_data$plantOccurs)
  
  
  # check relative size difference datasets
  rel = species_plot_number/size
  
  if(rel< 0.10){
    needed_size <- (species_plot_number/0.10)-species_plot_number
    data_to_sample <- subset_data[!subset_data$plantOccurs,]
    sample <- sample(data_to_sample$PlotObservationID, needed_size, replace=F)
    update <- data_to_sample[data_to_sample$PlotObservationID %in% sample,]
    subset_data <- rbind(subset_data[subset_data$plantOccurs,], update)
    size<- length(subset_data$plantOccurs)
  }
  
  if(rel> 0.30){
    needed_size <- round((species_plot_number/0.30)-size)
    data_to_sample <- tmpFullPlotData[!tmpFullPlotData$PlotObservationID %in% subset_data$PlotObservationID,]
    sample <- sample(data_to_sample$PlotObservationID, needed_size, replace=F)
    update <- data_to_sample[data_to_sample$PlotObservationID %in% sample,]
    subset_data <- rbind(subset_data, update)
    size<- length(subset_data$plantOccurs)
  }  
  
  if(size < 100){
    needed_size <- 100-size
    data_to_sample <- tmpFullPlotData[!tmpFullPlotData$PlotObservationID %in% subset_data$PlotObservationID,]
    sample <- sample(data_to_sample$PlotObservationID, needed_size, replace=F)
    update <- data_to_sample[data_to_sample$PlotObservationID %in% sample,]
    subset_data <- rbind(subset_data, update)
    size<- length(subset_data$plantOccurs)
  }
  
  
  # change name tmpFullPlotData
  tmpFullPlotData <- subset_data
  
  # If species is not present:
  # Make 1 total cover and 0 species cover
  tmpFullPlotData$totalCover[is.na(tmpFullPlotData$totalCover)] <- 100.0
  tmpFullPlotData$speciesCover[is.na(tmpFullPlotData$speciesCover)] <- 0.0
  
  # Exclude all plots which have 0 total cover
  tmpFullPlotData <- tmpFullPlotData[tmpFullPlotData$totalCover >0,] 
  
  # Do SR - 1 for all plots in which the species is present
  tmpFullPlotData$numberOfVascularPlantSpecies <- tmpFullPlotData$numberOfVascularPlantSpecies - tmpFullPlotData$plantOccurs
  
  correctArea= F
  # Correct area: if true we calculate the area that is reserved for the other plant species. 
  if(correctArea) {
    tmpFullPlotData$Area <- tmpFullPlotData$Area*(1.0 - tmpFullPlotData$speciesCover/tmpFullPlotData$totalCover)
  }
  
  # There is a numerical problem in the bam function if we provide empty plots...
  tmpFullPlotData <- tmpFullPlotData |> filter(numberOfVascularPlantSpecies > 0, Area > 0)
  
  # Make dataset a factor to add it as a random variable
  tmpFullPlotData$Dataset <- as.factor(tmpFullPlotData$Dataset)
  
  # Add relative species cover
  tmpFullPlotData <- tmpFullPlotData %>% 
    mutate(relSpeciesCover = 100 * speciesCover / totalCover)
  
  
  # Species occurrence as factor (not logical)
  tmpFullPlotData$plantOccursF <- factor(tmpFullPlotData$plantOccurs)
  
  # Base model
  model<-  numberOfVascularPlantSpecies ~ 
    s(log(Area),bs='tp') + 
    s(EIVEresM, bs = 'tp') +
    s(EIVEresN, bs = 'tp') +
    s(EIVEresR, bs = 'tp') +
    s(EIVEresL, bs = 'tp') +
    s(EIVEresT, bs = 'tp') +
    s(DistSeverity.sqrt, bs = 'tp') +
    s(Soil.Disturbance.sqrt, bs = 'tp')+
    s(Grazing.Pressure.sqrt, bs = 'tp')+
    s(Mowing.Frequency.sqrt, bs = 'tp') +
    s(Latitude, Longitude, bs = 'tp') +
    s(hfp, bs='tp')+
    s(chelsaP, bs='tp')+
    s(Dataset, bs = 're')+
    plantOccursF
  
  base <- bam(model,family = poisson(link = log), 
              data = tmpFullPlotData,  method = 'fREML',  discrete=TRUE, 
              nthreads=4)
  
  # Get summary from model
  baseModel <- summary(base)
  baseModel
  
  # Make a dataframe to add the residuals to
  res <- data.frame(taxa = character(), Estimate = numeric(), 
                    StdErr = numeric(), `zValue` = numeric(), 
                    pr = numeric(), numberOfPlots = integer(), 
                    Neophyte= character(),
                    RelDiff= numeric(),
                    size= numeric(),
                    species_plots= numeric(),
                    rel= numeric())
  
  res <- add_row(res, 
                 taxa = plantName, Estimate = baseModel$p.table[2,1], 
                 StdErr = baseModel$p.table[2,2], 
                 zValue = baseModel$p.table[2,3], 
                 pr = baseModel$p.table[2,4], 
                 numberOfPlots = length(plotsWherePlantOccurs$PlotObservationID), 
                 Neophyte= plantStatus,
                 RelDiff= 100* 
                   (mean(fitted(base)[tmpFullPlotData$plantOccursF == TRUE]) -
                      mean(fitted(base)[tmpFullPlotData$plantOccursF == FALSE])) /
                   mean(fitted(base)[tmpFullPlotData$plantOccursF == FALSE]),
                 size= size, 
                 species_plots = species_plot_number,
                 rel =species_plot_number/size)
  
  res
}
```

#### BRMS ####
# perform fancy test (has to be updated)
if(fancy){  
  library(brms)
  # bayesian model test with both mean and sigma depending on the group
  fit1 <- brm(bf(RelDiff ~ Neophyte, sigma ~ Neophyte, alpha ~Neophyte), data = Dataset, 
              family = "skew_normal", warmup= 1000, iter= 10000, chains= 4, cores = 12, control= list(adapt_delta=0.90))
  summary(fit1)
  
  # calculate the differences between the neophytes
  ggemmeans(fit1,"Neophyte",x.as.factor = TRUE)
  
  # test hypothesis (does not work)
  hypothesis(fit1, "Neophyteintra - Intercept=0")
  
  # also way to compare the species
  fit2.emm.a <- emmeans(fit1, "Neophyte", data=Dataset)
  pairs(fit2.emm.a, adjust="tukey")
  plot(fit2.emm.a, comparisons = F)
  
  # make 0 models (both with and without sigma)
  fit2 <- brm(RelDiff~1,data=Dataset)
  fit3 <- brm(bf(RelDiff~1, sigma ~Neophyte),data=Dataset)
  
  # calculate Bayes factor
  bayes_factor(fit1,fit3)
  
  #Check model
  plot(fit1)
  pp_check(fit1)
  summary(fit1)
}  


#### ART ####
result$eta <- with(result, `Sum Sq`/(`Sum Sq` + `Sum Sq.res`))
result
Sum = groupwiseMedian(RelDiff ~ Neophyte*class,
                      data=x,
                      bca=FALSE, percentile=TRUE)
Sum

pd = position_dodge(0.4)

ggplot(Sum,
       aes(x     = class,
           y     = Median,
           color = Neophyte)) +
  geom_point(shape  = 16,
             size   = 2,
             position = pd) +
  geom_errorbar(aes(ymin  =  Percentile.lower,
                    ymax  =  Percentile.upper),
                width =  0.2,
                size  =  0.7,
                position = pd) +
  theme_bw() +
  theme(axis.title   = element_text(face = "bold"),
        axis.text    = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0))



#### TRY  ####


# 2 TRY
## 2.1 Explore
```{r}
# Get all trait data
Traits <- rtry_import("../TRY/32370.txt", separator = "\t",encoding = "Latin-1", quote = "",showOverview = TRUE)

# How is data structured
head(Traits)
colnames(Traits)

# Which traits are in the dataset
Trait_variability<- rtry_explore(Traits, TraitID, TraitName)
Species_variability<- rtry_explore(Traits, AccSpeciesID, AccSpeciesName, TraitID, TraitName)
```


## 2.2 Select
```{r}
# select columns
workdata <- rtry_select_col(Traits, ObsDataID, ObservationID, AccSpeciesID, AccSpeciesName, 
                            ValueKindName, TraitID, TraitName, DataID, DataName, OriglName, 
                            OrigValueStr, OrigUnitStr, StdValue, UnitName, OrigObsDataID, 
                            ErrorRisk, Comment)

# explore sorted dataframe
workdata_explore_anc <- rtry_explore(workdata, DataID, DataName, TraitID, TraitName, sortBy = TraitID)

# Data ID
# 659, 60, 61    long, lat, alt
# 62, 80         MAT, MAP

# select which rows we want to keep
workdata <- rtry_select_row(workdata, TraitID > 0 | DataID %in% c(59, 60, 61, 62, 80, 413))
```


```{r}
# Save a version before deleting data
workdata_unexcluded <- workdata
```


## 2.3 Exclude
Select rows we want to work with
```{r}
# explore dataframe
tmp_unfiltered <- rtry_explore(workdata, DataID, DataName, TraitID, TraitName, sortBy = TraitID)

# exclude some data
workdata <- rtry_exclude(workdata, DataID %in% c(974, 1629, 1739, 3727, 3728, 8178, 9780:9782,8681,8682), baseOn = ObsDataID)

# select some rows
all<- rtry_select_row(workdata, DataID %in% c(4, 2568, 14, 15, 16, 30, 31, 258, 65, 100, 407, 
                                              485, 620, 19,20,448,504, 6575, 6463,6577,6579, 
                                              6581, 6589,6582, 6584,6598), baseOn=ObsDataID)

# explore newly made dataframe
tmp_unfiltered <- rtry_explore(all, DataID, DataName, TraitID, TraitName, sortBy = TraitID)
```


Remove values that are likely errors
```{r}
# error risk --> distance species/ genus mean --> 3 --> 3 sd away from mean
tmp_unfiltered <- rtry_explore(all, DataID, DataName, TraitID, TraitName, ErrorRisk, sortBy = ErrorRisk)
all <- rtry_exclude(all, ErrorRisk >= 5, baseOn = ObsDataID)

```


Duplicates
```{r}
# Based on OrigObsDataID but important to note that also data is removed that might not be duplicated if not all data was available from the beginning
workdata <- rtry_remove_dup(workdata)
```


## 2.4 Transform
```{r}
# from long to wide (more columns, less rows)

# only select rows with numeric values
num_traits <- rtry_select_row(workdata, complete.cases(TraitID) & complete.cases(StdValue))
# take the columns we want to use
num_traits <- rtry_select_col(num_traits, ObservationID, AccSpeciesID, AccSpeciesName, TraitID, TraitName, StdValue, UnitName)

# before transformation summarise
num_traits <- num_traits |> group_by(AccSpeciesName, TraitID, TraitName, UnitName) |> summarise(StdValue= mean(StdValue))

# transformation --> get names from trait names and use values as cell values (with mean)
num_traits_wider <-rtry_trans_wider(num_traits, names_from = c(TraitID, TraitName, UnitName), values_from = c(StdValue), values_fn = list(StdValue = mean))
```


## 2.5 Eva
```{r}
# load data
eva <- read_csv("fullPlotEva_cover_all_layer_cleaned.csv",show_col_types = FALSE)

# get species names of both datasets
species_try <- unique(num_traits_wider$AccSpeciesName)
species_eva<- unique(eva$species)
# count how many of our species are present in TRY
sum(species_eva %in% species_try)

# merge eva and trait data
test<- left_join(eva, num_traits_wider, by=c("species"="AccSpeciesName"))

# Count number of NA values per trait
na_counts <- (apply(test[, 21:48], 2, function(x) sum(!is.na(x))))
na_counts_df <- data.frame(variable = names(na_counts), value = na_counts)
na_counts_df

# Summarise data per plot
x<- test |> group_by(PlotObservationID) |> 
  summarise(cwm_SSD= mean(`4_Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density_g/cm3`, na.rm=T), 
            cwm_SLA= mean(`3115_Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded_mm2 mg-1`, na.rm=T), 
            n=n(`4_Stem specific density (SSD, stem dry mass per stem fresh volume) or wood density_g/cm3`))
```

#### PGLS ####
```{r}
library(nlme)
library(ape)
library(caper)

mat <- vcv(tree$scenario.3, corr=TRUE)

colnames(impact)[1] <- c("Names")
```


```{r}
data <- comparative.data(phy= tree$scenario.3, data= impact)


fit <- pgls(impact~ LMA, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=general, na.action= na.omit)
summary(fit)


```

```{r}
mat <- vcv(tree$scenario.3, corr=TRUE)

fit <- gls(RelDiff~ LMA, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=impact, na.action= na.omit)
summary(fit)

impact3 <- impact[!is.na(impact$RelDiff),]
impact3 <- impact3[!is.na(impact3$LMA),]
pglsModel <- gls(RelDiff ~ LMA , correlation = corPagel(1, phy = tree$scenario.3),  data = impact3, method = "ML")
summary(pglsModel)


gls(RelDiff ~ LMA, correlation = corBrownian(phy = tree$scenario.3, form= ~Names),data= impact, method="ML")
```

#### Significance test dominance ####

speciesDominance<- speciesDominance[speciesDominance$coverMean>0,]
speciesDominance <- speciesDominance[!speciesDominance$names %in% genus$original,]

# test the distributions
hist(log10(speciesDominance$coverMean))
hist(log10(speciesDominance$coverMean[speciesDominance$neophyte=="extra"]))
hist(log10(speciesDominance$coverMean[speciesDominance$neophyte=="native"]))
hist(log10(speciesDominance$coverMean[speciesDominance$neophyte=="intra"]))
#hist(log10(speciesDominance$coverMean[speciesDominance$neophyte=="native_intra"]))

# test the variance per group 
leveneTest(log10(speciesDominance$coverMean), speciesDominance$neophyte, location="mean")
# not normal residuals
shapiro.test(resid(aov(log10(speciesDominance$coverMean) ~ speciesDominance$neophyte))[1:5000])
hist(resid(aov(log10(speciesDominance$coverMean) ~ speciesDominance$neophyte)))
# test homoscedasticity
bartlett.test(log10(speciesDominance$coverMean), speciesDominance$neophyte)


```{r}
# non parametric

# ranked test
model<-(aov(rank(log10(coverMean)) ~ neophyte, data= speciesDominance))
summary(model)
TukeyHSD(model)
plot(TukeyHSD(model, conf.level=.95), las = 2)

# Welch Anova
model<-(oneway.test((log10(coverMean)) ~ neophyte, data= speciesDominance, var.equal = F))
model
tt1 <- posthoc_anova(log10(coverMean) ~ neophyte, data= speciesDominance, method="Tukey")
#comp<- agricolae::kruskal(log10(speciesDominance$coverMean),speciesDominance$neophyte, p.adj="BH")
tt1$groups

# kruskal wallis test
KW<-kruskal.test(log10(coverMean) ~ neophyte, data = speciesDominance)
KW$p.value
pairwise.wilcox.test(log10(speciesDominance$coverMean), speciesDominance$neophyte,
                     p.adjust.method = "BH")
dunnTest(log10(coverMean) ~ neophyte, data = speciesDominance,method = "bh")
```


# distribution
hist((log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="extra"])))
hist(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="native"]))
hist(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="intra"]))
# test the variance per group 
leveneTest(log10(speciesDominance$numberOfPlots), speciesDominance$neophyte, location="mean")
# not normal residuals
shapiro.test(resid(aov(log10(speciesDominance$numberOfPlots) ~ speciesDominance$neophyte))[1:5000])
hist(resid(aov(log10(speciesDominance$numberOfPlots) ~ speciesDominance$neophyte)))
# test homoscedasticity
bartlett.test(log10(speciesDominance$numberOfPlots), speciesDominance$neophyte)


# ranked test
model<-(aov(log10(numberOfPlots) ~ neophyte, data= speciesDominance))
summary(model)
TukeyHSD(model)
plot(TukeyHSD(model, conf.level=.95), las = 2)


# kruskal wallis test
kruskal.test(log10(numberOfPlots) ~ neophyte, data = speciesDominance)
pairwise.wilcox.test(log10(speciesDominance$numberOfPlots), speciesDominance$neophyte,
                     p.adjust.method = "BH")
dunnTest(log10(coverMean) ~ neophyte, data = speciesDominance,method = "bh")


##### Dominance plot ####
if(!native_intra_analysis){
  my_comparisons <- list( c("extra", "intra"), c("intra", "native"), c("extra", "native") )
  
  p<-ggplot(speciesDominance, mapping = aes(x= (neophyte), y = log10(coverMean), group= neophyte, colour = neophyte, fill= neophyte)) + 
    # create violin plot with x axis no title and y 
    geom_violin(alpha= 0.5, scale= "width")+
    geom_boxplot(width= 0.25, alpha=0.8, fill="white") +
    theme_pubr()+
    stat_compare_means(comparisons=my_comparisons, label= "p.signif", label.y = c(2.1, 2.3, 2.5))+
    ylab("Mean cover when present (%)") +
    scale_colour_manual(values=c("#1E88E5", "#FFC107", "#004D40"), 
                        breaks=c("native", "intra","extra"),
                        labels=c("native species", "intra European aliens", "extra European aliens")) +
    guides(color="none")+
    theme(legend.position = "none")+
    stat_summary(fun= "mean",
                 geom = "point", aes(group= neophyte), size=3)+
    scale_y_continuous(breaks = c(-1, 0, 1, 2), labels = c(0.1,1,10,100), limits=c(-1, 2.7))+
    scale_x_discrete(labels=  c("extra-European \n aliens",  "intra-European \n aliens", "native species")) +
    theme(axis.text.x = element_text(size=12), axis.title.y=element_text(size=12), axis.title.x = element_blank(), 
          axis.text.y = element_text(size=10)) +
    annotate("text", x = 1, y = mean(log10(speciesDominance$coverMean[speciesDominance$neophyte=="extra"])), label = "a", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8) +
    annotate("text", x = 2, y = mean(log10(speciesDominance$coverMean[speciesDominance$neophyte=="intra"])), label = "b", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x = 3, y = mean(log10(speciesDominance$coverMean[speciesDominance$neophyte=="native"])), label = "c", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x=1, y=2.0, label= paste("n=",sum(speciesDominance$neophyte=="extra")), vjust = 0, hjust = 0.5, size=3.5)+
    annotate("text", x=2, y=2.0, label= paste("n=",sum(speciesDominance$neophyte=="intra")), vjust = 0, hjust = 0.5, size=3.5)+
    annotate("text", x=3, y=2.0, label= paste("n=",sum(speciesDominance$neophyte=="native")), vjust = 0, hjust = 0.5, size=3.5)+
    labs(subtitle = substitute(paste("Kruskal-Wallis test ", italic("P < 0.0001"))))
  Dominance <- p
  
  #ggsave("Dominance New 30.jpeg", p,width = 5, height = 5)
}


if(native_intra_analysis){
  my_comparisons <- list( c("extra", "intra"), c("intra", "native"), c("extra", "native"), c("native", "native_intra"),
                          c("native_intra", "intra"), c("native_intra", "extra"))
  
  p<-ggplot(speciesDominance, mapping = aes(x= neophyte, y = log10(coverMean), group= neophyte,  fill= neophyte, colour=neophyte)) + 
    # create violin plot with x axis no title and y 
    geom_violin(alpha= 0.3, scale= "width")+
    geom_boxplot(width= 0.25, alpha=0.8, fill="white") +
    theme_pubr()+
    stat_compare_means(comparisons=my_comparisons, label= "p.signif", label.y = c(2.1, 2.3, 2.5, 2.7,2.9,3.1), size= 3)+
    ylab("Mean cover when present (%)") +
    scale_colour_manual(values=c( "#004D40","#FFC107", "#1E88E5", "darkgreen"), 
                        breaks=c("extra","intra", "native",  "native_intra")) +
    scale_fill_manual(values=c( "#004D40","#FFC107", "#1E88E5", "darkgreen"), 
                      breaks=c("extra","intra", "native",  "native_intra")) +
    guides(color="none")+
    theme(legend.position = "none")+
    stat_summary(fun= "mean",
                 geom = "point", aes(group= neophyte), size=3)+
    scale_y_continuous(breaks = c(-1, 0, 1, 2), labels = c(0.1,1,10,100), limits=c(-1, 3.3))+
    scale_x_discrete(labels=  c("extra-European \n neophytes",  "intra-European \n neophytes","native species", 
                                "native species \nalien elsewhere")) +
    theme(axis.text.x = element_text(size=10), axis.title.y=element_text(size=12), axis.title.x = element_blank(), 
          axis.text.y = element_text(size=10)) +
    annotate("text", x = 1, y = mean(log10(speciesDominance$coverMean[speciesDominance$neophyte=="extra"])), label = "a", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8) +
    annotate("text", x = 2, y = mean(log10(speciesDominance$coverMean[speciesDominance$neophyte=="intra"])), label = "b", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x = 3, y = mean(log10(speciesDominance$coverMean[speciesDominance$neophyte=="native"])), label = "d", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x = 4, y = mean(log10(speciesDominance$coverMean[speciesDominance$neophyte=="native_intra"])), label = "c", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x=1, y=2.0, label= paste("n=",sum(speciesDominance$neophyte=="extra")), vjust = 0, hjust = 0.5, size=3.5)+
    annotate("text", x=2, y=2.0, label= paste("n=",sum(speciesDominance$neophyte=="intra")), vjust = 0, hjust = 0.5, size=3.5)+
    annotate("text", x=3, y=2.0, label= paste("n=",sum(speciesDominance$neophyte=="native")), vjust = 0, hjust = 0.5, size=3.5)+
    annotate("text", x=4, y=2.0, label= paste("n=",sum(speciesDominance$neophyte=="native_intra")), vjust = 0, hjust = 0.5, size=3.5)+
    labs(subtitle = substitute(paste("Kruskal-Wallis test ", italic("P < 0.0001"))))
  p
  Dominance <- p
}
p

p<- ggplot(speciesDominance, mapping = aes(x= neophyte, y = log10(numberOfPlots), colour = neophyte, group=neophyte, fill= neophyte)) + 
  geom_violin(alpha= 0.5, scale= "width")+
  geom_boxplot(width= 0.25, alpha=0.8, fill="white") +
  theme_pubr()+
  stat_compare_means(comparisons=my_comparisons, label= "p.signif", label.y = c(5.5, 5.8, 6.1), size=3)+
  ylab(expression("log"[10]* " number of occurences")) + xlab(NULL)+
  scale_colour_manual(values=c("#1E88E5", "#FFC107", "#004D40"), 
                      name="Legend",
                      breaks=c("native", "intra","extra"),
                      labels=c("native species", "intra European aliens", "extra European aliens")) +
  guides(color="none")+
  stat_summary(fun= "mean",
               geom = "point", aes(group= neophyte), size=3)+
  guides(color="none")+
  theme(legend.position = "none")+
  scale_x_discrete(labels= c("extra-European \n aliens",  "intra-European \n aliens", "native species")) +
  theme(axis.text.x = element_text(size=12), axis.title.y=element_text(size=12) )  +
  theme(plot.subtitle=element_text(size=12)) + 
  annotate("text", x = 1, y = mean(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="extra"])), label = "b", size = 4, 
           vjust = -0.7, hjust = 0.5, alpha=0.8) +
  annotate("text", x = 2, y = mean(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="intra"])), label = "c", size = 4, 
           vjust = -0.7, hjust = 0.5, alpha=0.8)+
  annotate("text", x = 3, y = mean(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="native"])), label = "a", size = 4, 
           vjust = -0.7, hjust = 0.5, alpha=0.8)+
  annotate("text", x=1, y=5.2, label= paste("n=",sum(speciesDominance$neophyte=="extra")), vjust = -0.5, hjust = 0.5, size=3.5)+
  annotate("text", x=2, y=5.2, label= paste("n=",sum(speciesDominance$neophyte=="intra")), vjust = -0.5, hjust = 0.5, size=3.5)+
  annotate("text", x=3, y=5.2, label= paste("n=",sum(speciesDominance$neophyte=="native")), vjust = -0.5, hjust = 0.5, size=3.5)+
  labs(subtitle = substitute(paste("Kruskal-Wallis test ", italic("P < 0.0001"))))
Frequency <- p
p
#ggsave("Frequency New 10.jpeg", p,width = 6, height = 4)
} else {
  p<- ggplot(speciesDominance, mapping = aes(x= neophyte, y = log10(numberOfPlots), colour = neophyte, group=neophyte, fill= neophyte)) + 
    geom_violin(alpha= 0.3, scale= "width")+
    geom_boxplot(width= 0.25, alpha=0.8, fill="white") +
    theme_pubr()+
    stat_compare_means(comparisons=my_comparisons, label= "p.signif", label.y = c(5.5, 5.8, 6.1, 6.4,6.7,7), size=3)+
    ylab(expression("log"[10]* " number of occurences")) + xlab(NULL)+
    scale_colour_manual(values=c("#1E88E5", "#FFC107", "#004D40", "darkgreen"), 
                        name="Legend",
                        breaks=c("native", "intra","extra", "native_intra"),
                        labels=c("native species", "intra European aliens", "extra European aliens", "native species alien elsewhere")) +
    scale_fill_manual(values=c( "#004D40","#FFC107", "#1E88E5", "darkgreen"), 
                      breaks=c("extra","intra", "native",  "native_intra")) +
    guides(color="none")+
    stat_summary(fun= "mean",
                 geom = "point", aes(group= neophyte), size=3)+
    guides(color="none")+
    theme(legend.position = "none")+
    scale_x_discrete(labels= c("extra-European \n aliens",  "intra-European \n aliens", "native species",
                               "native species \nalien elsewhere")) +
    theme(axis.text.x = element_text(size=10), axis.title.y=element_text(size=12) )  +
    theme(plot.subtitle=element_text(size=12)) + 
    annotate("text", x = 1, y = mean(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="extra"])), label = "c", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8) +
    annotate("text", x = 2, y = mean(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="intra"])), label = "d", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x = 3, y = mean(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="native"])), label = "b", size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x = 4, y = mean(log10(speciesDominance$numberOfPlots[speciesDominance$neophyte=="native_intra"])), label = "a",size = 4, 
             vjust = -0.7, hjust = 0.5, alpha=0.8)+
    annotate("text", x=1, y=5.2, label= paste("n=",sum(speciesDominance$neophyte=="extra")), vjust = -0.5, hjust = 0.5, size=3.5)+
    annotate("text", x=2, y=5.2, label= paste("n=",sum(speciesDominance$neophyte=="intra")), vjust = -0.5, hjust = 0.5, size=3.5)+
    annotate("text", x=3, y=5.2, label= paste("n=",sum(speciesDominance$neophyte=="native")), vjust = -0.5, hjust = 0.5, size=3.5)+
    annotate("text", x=4, y=5.2, label= paste("n=",sum(speciesDominance$neophyte=="native_intra")), vjust = -0.5, hjust = 0.5, size=3.5)+
    labs(subtitle = substitute(paste("Kruskal-Wallis test ", italic("P < 0.0001"))))
  Frequency <- p
  p  

  
  #### BRTs ####
  
  ## 2.5 BRTs
  Try boosted regression tree
  ```{r}
  library(gbm)
  library(dismo)
  
  fullPlotData2 <- as.data.frame(fullPlotData)
  BRTmodel1 <- gbm.step(data=fullPlotData2,
                        gbm.x = c(4:12, 27:28), 
                        gbm.y = 2, 
                        family = "poisson",
                        tree.complexity = 4,
                        learning.rate = 0.001,
                        bag.fraction = 0.7)
  
  BRTmodel1$contributions
  gbm.plot(BRTmodel1, n.plots=11, write.title = F,
           common.scale = F, y.label="Fitted function",
           show.contrib = TRUE)
  
  BRTmodel1$contributions
  gbm.plot(BRTmodel1, n.plots=11, write.title = F,
           common.scale = F, y.label="Fitted function",
           show.contrib = TRUE)
  ```
#### GLMM ####

  
  
  ## 2.4 GLMM
  ```{r}
  # run model
  
  sqrt =F
  if(sqrt){
    MDL <- glmmTMB(ENS0 ~ EIVEresN+ hfp+ cover_rel_intra + cover_rel_extra+ EIVEresM+ EIVEresR+
                     EIVEresL+ EIVEresT+ DistSeverity + I(DistSeverity ^2)+ log(Area)+
                     Soil.Disturbance + I(Soil.Disturbance ^2)+ 
                     Grazing.Pressure +I(Grazing.Pressure ^2)+ Mowing.Frequency +
                     I(Mowing.Frequency ^2)+ chelsaP+Longitude+Latitude+ (1| Dataset), 
                   data=fullPlotData, family= poisson(link=log))
  } else {
    MDL <- glmmTMB(n_native ~ EIVEresN+ hfp+ n_intra+ n_extra + EIVEresM+ EIVEresR+
                     EIVEresL+ EIVEresT+ DistSeverity.sqrt + I(DistSeverity.sqrt ^2)+ log(Area)+
                     Soil.Disturbance.sqrt + I(Soil.Disturbance.sqrt ^2)+ 
                     Grazing.Pressure.sqrt +I(Grazing.Pressure.sqrt ^2)+ Mowing.Frequency.sqrt +
                     I(Mowing.Frequency.sqrt ^2)+ chelsaP+Longitude+Latitude+ (1| Dataset), 
                   data=fullPlotData, family= poisson(link=log))
  }
  
  
  # the sqrt model works better
  sum<- summary(MDL)
  sum
  piecewiseSEM::rsquared(MDL)
  ```
  
  
  
  
  
  ```{r}
  # plot using visreg
  z<- visreg::visreg(MDL, xvar="EIVEresN", gg = TRUE, partial=FALSE, rug = F,  scale="response", plot=F)
  eive_n_coef <- sum$coefficients$cond[2,1]/ z$fit$visregFit[1]
  
  z2 <- visreg::visreg(MDL, xvar="n_extra", gg = TRUE, partial=FALSE, rug = F,  scale="response", plot=F)
  extra_coef <- sum$coefficients$cond[5,1]/ z2$fit$visregFit[1]
  
  z3 <-visreg::visreg(MDL, xvar="n_intra", gg = TRUE, partial=FALSE, rug = F,  scale="response", plot=F)
  intra_coef <- sum$coefficients$cond[4,1]/ z3$fit$visregFit[1]
  ```
  
  
  
  Check
  ```{r}
  # check data
  library(DHARMa)
  
  # create simulation data
  simulationOutput <- simulateResiduals(fittedModel = MDL, plot = F)
  # look at simulation 
  # QQ plot --> line normal (is it correct distribution and is there no overdispersion)
  # residual plot --> if correct line at 0.5
  plot(simulationOutput)
  
  # test dispersion
  testDispersion(simulationOutput)
  plotResiduals(simulationOutput)
  
  # are there more outliers than expected by accident
  testOutliers(simulationOutput)
  
  # compare location with expected
  testQuantiles(simulationOutput)
  
  # more zeroes than expected
  testZeroInflation(simulationOutput)
  
  # plot per random effect group
  simulationOutput2 <- recalculateResiduals(simulationOutput, group = fullPlotData$Dataset)
  plot(simulationOutput2, quantreg = FALSE) 
  ```
  
  ```{r}
  library(xgboost)
  train <- fullPlotData2[runif(length(fullPlotData2$PlotObservationID)) > 0.30,]
  test <- fullPlotData2[!fullPlotData2$PlotObservationID %in% train$PlotObservationID,]
  
  bst <- xgboost(data = as.matrix(train[, c(4:12)]), label = train$ENS0, max.depth = 4, eta = 1, nthread = 1, nrounds = 20, objective = "count:poisson")
  
  pred <- predict(bst, as.matrix(test[,c(4:12)]))
  err <- mean(as.numeric(pred)- test$ENS0)
  err
  max(pred)
  mean(pred)
  max(train$ENS0)
  mean(train$ENS0)
  
  ```
  
  
  
##### EIVE AND DIV #####  
  ## 2.1 Eive
  Reading in the EIVE data and explore data by boxplotting the moisture (M), nitrogen (N), pH (R), light (L) and temperature (T) variables.
  ```{r, message=F, fig.show = "hide"}
  # 08 contains the indicator values for position and width
  eive <- read_delim("../Extra data/ENV/EIVE_Paper_1.0_SM_08.csv", ",")
  # data contains n --> number of sources, nw--> niche width and ''--> value
  
  par(mfrow = c(1,1))
  boxplot(eive$`EIVEres-M`, xlab="M") # outliers above
  boxplot(eive$`EIVEres-N`, xlab="N") # little outliers
  boxplot(eive$`EIVEres-R`, xlab="R") # outliers below
  boxplot(eive$`EIVEres-L`, xlab="L") # outliers below
  boxplot(eive$`EIVEres-T`, xlab="T") # little outliers
  par(mfrow=c(1,1))
  # The EIVE values are OK balanced with some outliers. No extreme skewness to be expected for the average indicator values at the plot level.
  ```
  
  
  Extend eva with columns from EIVE.
  ```{r, fig.show = "hide"}
  # match gives position of species in eva in eive
  
  # First we match our own 'species column', which has the Euro Med data added with powo and some data from Irena
  # Later 'the matched concept' was added
  
  eva$eive_name <- eive$TaxonConcept[match(eva$name, eive$TaxonConcept)]
  length(unique(eva$eive_name[!is.na(eva$eive_name)]))
  eva$eive_name[is.na(eva$eive_name)] <- eive$TaxonConcept[match(gsub(" aggr\\.", "", eva$name[is.na(eva$eive_name)]),
                                                                 eive$TaxonConcept)]
  length(unique(eva$eive_name[!is.na(eva$eive_name)]))
  eva$eive_name[is.na(eva$eive_name)] <- eive$TaxonConcept[match(eva$species[is.na(eva$eive_name)] , eive$TaxonConcept)]
  length(unique(eva$eive_name[!is.na(eva$eive_name)]))
  eva$eive_name[is.na(eva$eive_name)] <- eive$TaxonConcept[match(eva$irena[is.na(eva$eive_name)] , eive$TaxonConcept)]
  length(unique(eva$eive_name[!is.na(eva$eive_name)]))
  eva$eive_name[is.na(eva$eive_name)] <- eive$TaxonConcept[match(eva$`Matched concept`[is.na(eva$eive_name)],
                                                                 eive$TaxonConcept)]
  length(unique(eva$eive_name[!is.na(eva$eive_name)]))
  eva$eive_name[is.na(eva$eive_name)] <- eive$TaxonConcept[match(eva$`Turboveg2 concept`[is.na(eva$eive_name)],
                                                                 eive$TaxonConcept)]
  length(unique(eva$eive_name[!is.na(eva$eive_name)]))
  
  # positions
  eva$EIVEresM <- eive$`EIVEres-M`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEresN <- eive$`EIVEres-N`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEresR <- eive$`EIVEres-R`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEresL <- eive$`EIVEres-L`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEresT <- eive$`EIVEres-T`[match(eva$eive_name, eive$TaxonConcept)]
  
  # niche widths
  eva$EIVEnwM <- eive$`EIVEres-M.nw3`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEnwN <- eive$`EIVEres-N.nw3`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEnwR <- eive$`EIVEres-R.nw3`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEnwL <- eive$`EIVEres-L.nw3`[match(eva$eive_name, eive$TaxonConcept)]
  eva$EIVEnwT <- eive$`EIVEres-T.nw3`[match(eva$eive_name, eive$TaxonConcept)]
  ```
  
  
  
  ```{r, fig.show = "hide"}
  # look at total number of species in eive
  length(unique(eive$TaxonConcept[match(gsub(" aggr\\.", "", eva$name), eive$TaxonConcept)]))
  length(unique(eive$TaxonConcept[match(eva$species, eive$TaxonConcept)]))
  length(unique(eive$TaxonConcept[match(eva$irena, eive$TaxonConcept)]))
  length(unique(eive$TaxonConcept[match(eva$`Matched concept`,eive$TaxonConcept)]))
  length(unique(eive$TaxonConcept[match(eva$`Turboveg2 concept`,eive$TaxonConcept)]))
  
  
  # check total amount of observations
  sum(!is.na(eva$species[match(gsub(" aggr\\.", "", eva$name), eive$TaxonConcept)]))
  sum(!is.na(eva$species[match(eva$species, eive$TaxonConcept)]))
  sum(!is.na(eva$species[match(eva$irena, eive$TaxonConcept)]))
  sum(!is.na(eva$species[match(eva$`Matched concept`, eive$TaxonConcept)]))
  sum(!is.na(eva$species[match(eva$`Turboveg2 concept`, eive$TaxonConcept)]))
  
  # no names
  no_eive <-  eva[is.na(eva$eive_name),]
  no_eive <- unique(no_eive[, 2:6])
  
  # check against the length of eva names
  eva_names <- unique(eva[,2:6])
  
  percentOfSpeciesEvaInEive <- length(unique(eva$name[!is.na(eva$eive_name)]))/ length(unique(eva$name))
  
  
  # Look at number of observations not in EIVE
  numberOfPlantsInEvaFoundInEive <- sum(!(is.na(eva$EIVEresM) | 
                                            is.na(eva$EIVEresN) | 
                                            is.na(eva$EIVEresR) | 
                                            is.na(eva$EIVEresL) | 
                                            is.na(eva$EIVEresT)))
  percentEvaInEive<- numberOfPlantsInEvaFoundInEive/length(eva$PlotObservationID)
  percentEvaInEive
  percentOfSpeciesEvaInEive
  ```
  
  
  
  
  
  `r numberOfPlantsInEvaFoundInEive` of the `r length(eva$PlotObservationID)` plant observations in Eva have indicator values in EIVE.
  
  
  Look at unique species to see whether neophytes are not with less relative representatives
  This has to be performed with old data, since new data from Irena (January 2024) is only on the full plot data (which is made later in this file, after removing some plots etc)
  ```{r}
  # Load file
  taxonNameCorrections <- read_csv("../Extra data/Species names/Irena_taxa_check_Christian.csv", show_col_types = FALSE)
  
  # downsize to ease calculations
  test<- eva[, c("species", "irena","EIVEresM","EIVEresN","EIVEresR","EIVEresL","EIVEresT")]
  # remove duplicates species
  test<- test[!duplicated(test$species),]
  # join with taxonNameCorrections
  test<- left_join(test, taxonNameCorrections, by=c("irena"="species"))
  
  test$statusEurope<- as.character(test$statusEurope)
  # percent of natives without data
  length(test$species[(test$statusEurope=="native") & (is.na(test$EIVEresL)|                                                is.na(test$EIVEresM)|is.na(test$EIVEresN)|is.na(test$EIVEresR)|is.na(test$EIVEresT))])/length(test$species[test$statusEurope=="native"])
  
  # percent of neophytes without data
  length(test$species[(test$statusEurope=="neo") & (is.na(test$EIVEresL)| is.na(test$EIVEresM)|is.na(test$EIVEresN)|is.na(test$EIVEresR)|is.na(test$EIVEresT))])/length(test$species[test$statusEurope=="neo"])
  
  # percent of archeotypes without data
  length(test$species[(test$statusEurope=="arch") & (is.na(test$EIVEresL)| is.na(test$EIVEresM)|is.na(test$EIVEresN)|is.na(test$EIVEresR)|is.na(test$EIVEresT))])/length(test$species[test$statusEurope=="arch"])
  ```
  
  
  
  
  
  
  ## 2.2 Div
  Reading in the disturbance indicator values en plotting severity and frequency.
  ```{r, messagge=F, warning=F, message=F, fig.show = "hide"}
  div <- read_delim("../Extra data/DIV/disturbance_indicator_values.csv", ",")
  
  par(mfrow = c(1,1))
  # Relatively balanced
  boxplot(div$Disturbance.Severity, xlab="Severity")
  boxplot(div$Disturbance.Frequency, xlab="Frequency")
  
  # more skewed (check with hist)
  boxplot(div$Mowing.Frequency, xlab="Mowing fr")
  hist(div$Mowing.Frequency, xlab="Mowing fr")
  
  # check others
  boxplot(div$Grazing.Pressure, xlab="Grazing Pr")
  boxplot(div$Soil.Disturbance, xlab="Soil dist")
  
  # oke balanced still
  hist(div$Grazing.Pressure, xlab="Grazing Pr")
  ```
  
  
  Extend eva with columns from disturbance indicator values:
    ```{r}
  # merge names
  eva$div_name <- div$species[match(eva$name, div$species)]
  length(unique(eva$div_name[!is.na(eva$div_name)]))
  eva$div_name[is.na(eva$div_name)] <- div$species[match(gsub(" aggr\\.", "",eva$name[is.na(eva$div_name)]), div$species)]
  length(unique(eva$div_name[!is.na(eva$div_name)]))
  eva$div_name[is.na(eva$div_name)] <- div$species[match(eva$species[is.na(eva$div_name)], div$species)]
  length(unique(eva$div_name[!is.na(eva$div_name)]))
  eva$div_name[is.na(eva$div_name)] <- div$species[match(eva$irena[is.na(eva$div_name)], div$species)]
  length(unique(eva$div_name[!is.na(eva$div_name)]))
  eva$div_name[is.na(eva$div_name)] <- div$species[match(eva$`Matched concept`[is.na(eva$div_name)], div$species)]
  length(unique(eva$div_name[!is.na(eva$div_name)]))
  eva$div_name[is.na(eva$div_name)] <- div$species[match(eva$`Turboveg2 concept`[is.na(eva$div_name)], div$species)]
  length(unique(eva$div_name[!is.na(eva$div_name)]))
  
  # expand first with our species classification
  eva$Disturbance.Severity <- div$Disturbance.Severity[match(eva$div_name, div$species)]
  eva$Disturbance.Severity.herblayer <- div$Disturbance.Severity.herblayer[match(eva$div_name, div$species)]
  eva$Disturbance.Frequency <- div$Disturbance.Frequency[match(eva$div_name, div$species)]
  eva$Disturbance.Frequency.herblayer <- div$Disturbance.Frequency.herblayer[match(eva$div_name, div$species)]
  eva$Grazing.Pressure <- div$Grazing.Pressure[match(eva$div_name, div$species)]
  eva$Mowing.Frequency <- div$Mowing.Frequency[match(eva$div_name, div$species)]
  eva$Soil.Disturbance <- div$Soil.Disturbance[match(eva$div_name, div$species)]
  ```
  
  
  ```{r}
  # look at total number of species in eive
  length(unique(div$species[match(eva$name, div$species)]))
  length(unique(div$species[match(eva$species, div$species)]))
  length(unique(div$species[match(eva$irena, div$species)]))
  length(unique(div$species[match(eva$`Matched concept`,div$species)]))
  length(unique(div$species[match(eva$`Turboveg2 concept`,div$species)]))
  
  # check total amount of observations
  sum(!is.na(eva$species[match(eva$name, div$species)]))
  sum(!is.na(eva$species[match(eva$species, div$species)]))
  sum(!is.na(eva$species[match(eva$irena, div$species)]))
  sum(!is.na(eva$species[match(eva$`Matched concept`,div$species)]))
  sum(!is.na(eva$species[match(eva$`Turboveg2 concept`,div$species)]))
  
  # no names
  no_div <-  eva[is.na(eva$div_name),]
  no_div <- unique(no_div[, 2:4])
  
  # check against the length of eva names
  percentOfSpeciesEvaInDiv <- length(unique(eva$name[!is.na(eva$div_name)]))/ length(unique(eva$name))
  
  # and amount of observations
  numberOfPlantsInEvaFoundInDiv <- sum(!is.na(eva$Disturbance.Severity))
  percentEvaInDiv<- numberOfPlantsInEvaFoundInDiv/length(eva$Disturbance.Severity)
  percentEvaInDiv
  percentOfSpeciesEvaInEive
  ```
  
  `r numberOfPlantsInEvaFoundInDiv` of the `r length(eva$PlotObservationID)` plant observations in Eva have disturbance indicator value.
  
  Look at unique species to see whether neophytes are not with less relative representatives
  ```{r}
  # reduce eva
  test<- eva[, c("species","irena", "Disturbance.Severity")]
  # only unique species
  test<- test[!duplicated(test$species),]
  test<- left_join(test, taxonNameCorrections, by=c("irena"="species"))
  
  # make status a character
  test$statusEurope<- as.character(test$statusEurope)
  # percent of natives without data
  length(test$species[(test$statusEurope=="native") & (is.na(test$Disturbance.Severity))])/length(test$species[test$statusEurope=="native"])
  
  # percent of neophytes without data
  length(test$species[(test$statusEurope=="neo") & (is.na(test$Disturbance.Severity))])/length(test$species[test$statusEurope=="neo"])
  
  # percent of archeotypes without data
  length(test$species[(test$statusEurope=="arch") & (is.na(test$Disturbance.Severity))])/length(test$species[test$statusEurope=="arch"])
  ```
  
  
  
  ## 2.3 All
  ```{r}
  # create eva_names file again
  # vector with eva names
  eva_names<- as.data.frame(unique(eva$`Matched concept`))
  colnames(eva_names) <- "Matched concept"
  eva_names$name <- eva$name[match( eva_names$`Matched concept`, eva$`Matched concept`)]
  eva_names$species <- eva$species[match( eva_names$`Matched concept`, eva$`Matched concept`)]
  eva_names$eive <- eva$eive_name[match( eva_names$`Matched concept`, eva$`Matched concept`)]
  eva_names$div <- eva$div_name[match( eva_names$`Matched concept`, eva$`Matched concept`)]
  
  dup <- eva_names[duplicated(eva_names$name) | duplicated(eva_names$name, fromLast = T),]
  dup_known_eive <- dup[!is.na(dup$eive),]
  dup_known_div <- dup[!is.na(dup$div),]
  
  dup_unknown_eive <- dup[is.na(dup$eive),]
  dup_unknown_div <- dup[is.na(dup$div),]
  
  dup_unknown_eive <- dup_unknown_eive[dup_unknown_eive$name %in% dup_known_eive$name, ]
  dup_known_eive <- dup_known_eive[dup_known_eive$name %in% dup_unknown_eive$name,]
  
  dup_unknown_eive <- left_join(dup_unknown_eive[, -c(4)], dup_known_eive[, c(2,4)], by= c("name"="name"))
  ```
  
  
  
  EIVEresM.nw = NWwmean(EIVEresM, EIVEnwM), 
  EIVEresN.nw = NWwmean(EIVEresN, EIVEnwN), 
  EIVEresR.nw = NWwmean(EIVEresR, EIVEnwR), 
  EIVEresL.nw = NWwmean(EIVEresL, EIVEnwL),
  EIVEresT.nw = NWwmean(EIVEresT, EIVEnwT), 
  
  EIVEresM.cnw = cNWwmean(EIVEresM, `Cover %`, EIVEnwM), 
  EIVEresN.cnw = cNWwmean(EIVEresN, `Cover %`, EIVEnwN), 
  EIVEresR.cnw = cNWwmean(EIVEresR, `Cover %`, EIVEnwR), 
  EIVEresL.cnw = cNWwmean(EIVEresL, `Cover %`, EIVEnwL),
  EIVEresT.cnw = cNWwmean(EIVEresT, `Cover %`, EIVEnwT), 
  
  
  
  # 4 PLOT
  ## 4.1 SR
  ```{r}
  plot=F
  if(plot){
    par(mfrow=c(1,1))
    boxplot(fullPlotData$numberOfVascularPlantSpecies)
    hist(fullPlotData$numberOfVascularPlantSpecies)
    # As expected the count data on species richness looks quite skewed. Hence, we will need to use a poisson glm/ gam
  }
  ```
  
  
  ## 4.2 IVs
  ```{r}
  plot=F
  if(plot){
    # Square root
    boxplot(fullPlotData$DistSeverity.sqrt)
    boxplot(fullPlotData$DistFrequency.sqrt)
    boxplot(fullPlotData$Grazing.Pressure.sqrt)
    boxplot(fullPlotData$Mowing.Frequency.sqrt)
    boxplot(fullPlotData$Soil.Disturbance.sqrt)
    
    hist(fullPlotData$DistSeverity.sqrt)
    hist(fullPlotData$DistFrequency.sqrt)
    hist(fullPlotData$Grazing.Pressure.sqrt)
    hist(fullPlotData$Mowing.Frequency.sqrt)
    hist(fullPlotData$Soil.Disturbance.sqrt)
    
    boxplot(fullPlotData$EIVEresM.sqrt)
    boxplot(fullPlotData$EIVEresN.sqrt)
    boxplot(fullPlotData$EIVEresR.sqrt)
    boxplot(fullPlotData$EIVEresL.sqrt)
    boxplot(fullPlotData$EIVEresT.sqrt)
    
    # mean
    boxplot(fullPlotData$Grazing.Pressure)
    boxplot(fullPlotData$Mowing.Frequency)
    boxplot(fullPlotData$Soil.Disturbance)
    boxplot(fullPlotData$DistSeverity)
    boxplot(fullPlotData$DistFrequency)
    
    hist(fullPlotData$Grazing.Pressure)
    hist(fullPlotData$Mowing.Frequency)
    hist(fullPlotData$Soil.Disturbance)
    hist(fullPlotData$DistSeverity)
    hist(fullPlotData$DistFrequency)
    
    
    boxplot(fullPlotData$EIVEresM) # Looks good
    boxplot(fullPlotData$EIVEresN) # Looks good
    boxplot(fullPlotData$EIVEresR) # Looks good
    boxplot(fullPlotData$EIVEresL) # Looks good
    boxplot(fullPlotData$EIVEresT) # Looks good
    
    #Log mean
    boxplot(fullPlotData$logMowingFrequency)
    boxplot(fullPlotData$logGrazing.Pressure)
    boxplot(fullPlotData$logSoil.Disturbance)
    boxplot(fullPlotData$logDistSeverity)
    boxplot(fullPlotData$logDistFrequency)
    
    hist(fullPlotData$logMowingFrequency)
    hist(fullPlotData$logGrazing.Pressure)
    hist(fullPlotData$logSoil.Disturbance)
    hist(fullPlotData$logDistSeverity)
    hist(fullPlotData$logDistFrequency)
  }
  ```
  
  The distribution of the EIVE predictors look symmetrical. However, the disturbance indicators are heavily skewed. Hence, they need to be transformed. Because there are zeros in the data we choose a root transoformation. Trying different powers we found that the 4th root worked best:
    
    ## 4.3 Transform
    ```{r}
  plot=F
  if(plot){
    # originally 1/4
    fullPlotData$transformedDisturbanceSeverity.sqrt <- fullPlotData$DistSeverity.sqrt^(1/2)
    fullPlotData$transformedDisturbanceFrequency.sqrt <- fullPlotData$DistFrequency.sqrt^(1/2)
    fullPlotData$transformedDisturbanceSeverity <- fullPlotData$DistSeverity^(1/2)
    fullPlotData$transformedDisturbanceFrequency <- fullPlotData$DistFrequency^(1/2)
    
    hist(fullPlotData$transformedDisturbanceSeverity.sqrt)
    hist(fullPlotData$DistSeverity.sqrt) # imo even better distributed...
    hist(fullPlotData$transformedDisturbanceFrequency.sqrt)
    hist(fullPlotData$DistFrequency.sqrt)
    
    hist(fullPlotData$transformedDisturbanceSeverity)
    hist(fullPlotData$DistSeverity) # imo even better distributed...
    hist(fullPlotData$transformedDisturbanceFrequency)
    hist(fullPlotData$DistFrequency)
  }
  ```
  
  
  # 5 BASE MODEL
  Note that this is not the final Base Model, rather we will use this to check the assumptions made.
  ```{r}
  model=F
  if(model){
    glm0 <- glm(numberOfVascularPlantSpecies ~ 
                  log(Area) + 
                  EIVEresM + I(EIVEresM^2) + 
                  EIVEresN + I(EIVEresN^2) + 
                  EIVEresR + I(EIVEresR^2) + 
                  EIVEresL + I(EIVEresL^2) + 
                  EIVEresT + I(EIVEresT^2) +
                  DistSeverity + I(DistSeverity^2) +
                  DistFrequency + I(DistFrequency^2)
                , family=poisson, fullPlotData)
    
    summary(glm0)
    plot(glm0)
  }
  ```
  
  We can see that the residuals of this model are in a reasonable range and hence the model fits reasonably well to the data.
  ```{r}
  model=F
  if(model){
    dataFromGermany <- fullPlotData[fullPlotData$Region == "Germany",]
    dataFromGermany$Dataset <- as.factor(as.character(dataFromGermany$Dataset))
    par(mfrow = c(1,1))
    boxplot(numberOfVascularPlantSpecies ~ Dataset, data = dataFromGermany)
    boxplot(log(fullPlotData$Area))
  }
  ```
  
  Looking at different datasets from a similar region exhibits quite large differences in terms of species richness. It will most likely make sense to include Dataset as a random factor. The log transformed size of the plots is very symmetrical. This looks good.
  
  
  ##### Traits #####
  
  test <- hier[hier$Genus %in% duplicates$Genus,]
  
  test <- test |> group_by(Genus, Family) |> summarise(n=n())
  test <- test[order(test$n, decreasing=T),]
  test <- test[!duplicated(test$Genus),]
  
  
  ##### Traits 2 #####
  
  ## 4.5 Effect individual
  ```{r}
  # Make Neophyte a factor for in the models
  general$Neophyte <- as.factor(general$Neophyte)
  general$Neophyte <- relevel(general$Neophyte, ref="native")
  
  general[, c(3:8)] <- log(general[,c(3:8)])
  general[, c(3:8)] <- scale(general[, c(3:8)])
  
  # Model per FT individually, as limited species (+- 880) have values for all functional traits
  MDL_H<- lmer(impact ~ -1+Neophyte + H : Neophyte + (1|Family/Genus), general)
  summary(MDL_H)
  MDL_LMA<- lmer(impact ~ -1+Neophyte +LMA : Neophyte + (1|Family/Genus) , general)
  summary(MDL_LMA)
  MDL_SM<- lmer(impact ~ -1+Neophyte +SM : Neophyte + (1|Family/Genus), general)
  summary(MDL_SM)
  MDL_SSD<- lmer(impact ~-1+Neophyte + SSD : Neophyte + (1|Family/Genus), general)
  summary(MDL_SSD)
  MDL_N<- lmer(impact ~-1+Neophyte + N : Neophyte + (1|Family/Genus), general)
  summary(MDL_N)
  MDL_LA<- lmer(impact ~ -1+Neophyte +LA : Neophyte + (1|Family/Genus), general)
  summary(MDL_LA)
  
  
  # Plot
  visreg::visreg(MDL_SSD,"SSD", by="Neophyte")
  visreg::visreg(MDL_LMA,"LMA", by="Neophyte")
  visreg::visreg(MDL_LA,"LA", by="Neophyte")
  visreg::visreg(MDL_H,"H", by="Neophyte")
  visreg::visreg(MDL_N,"N", by="Neophyte")
  visreg::visreg(MDL_SM,"SM", by="Neophyte")
  ```
  
  
  
  And EIVE
  ```{r}
  # Make Neophyte a factor for in the models
  general$Neophyte <- as.factor(general$Neophyte)
  general$Neophyte <- relevel(general$Neophyte, ref="native")
  
  # scale EIVE?
  #general[, c(11:15)] <- scale(general[,c(11:15)])
  
  # Model per FT individually, as limited species (+- 880) have values for all functional traits
  MDL_T<- lmer(impact ~ -1+Neophyte + EIVEresT : Neophyte + (1|Family/Genus), general)
  summary(MDL_H)
  MDL_R<- lmer(impact ~ -1+Neophyte +EIVEresR : Neophyte + (1|Family/Genus), general)
  summary(MDL_LMA)
  MDL_N<- lmer(impact ~ -1+Neophyte +EIVEresN : Neophyte + (1|Family/Genus), general)
  summary(MDL_SM)
  MDL_L<- lmer(impact ~-1+Neophyte + EIVEresL : Neophyte + (1|Family/Genus), general)
  summary(MDL_SSD)
  MDL_M<- lmer(impact ~-1+Neophyte + EIVEresM : Neophyte + (1|Family/Genus), general)
  summary(MDL_N)
  
  
  
  # Plot
  visreg::visreg(MDL_M,"EIVEresM", by="Neophyte")
  visreg::visreg(MDL_L,"EIVEresL", by="Neophyte")
  visreg::visreg(MDL_N,"EIVEresN", by="Neophyte")
  visreg::visreg(MDL_R,"EIVEresR", by="Neophyte")
  visreg::visreg(MDL_T,"EIVEresT", by="Neophyte")
  ```
  
  
  And DIV
  ```{r}
  # Make Neophyte a factor for in the models
  general$Neophyte <- as.factor(general$Neophyte)
  general$Neophyte <- relevel(general$Neophyte, ref="native")
  general$woody <- as.factor(general$woody)
  general$woody <- relevel(general$woody,ref="woody")
  
  # scale EIVE?
  #general[, c(11:15)] <- scale(general[,c(11:15)])
  
  # Model per FT individually, as limited species (+- 880) have values for all functional traits
  MDL_sev<- lmer(impact ~ -1+Neophyte + Sev : Neophyte + (1|Family/Genus), general)
  summary(MDL_H)
  MDL_sev_herb<- lmer(impact ~ -1+Neophyte +Sev_herb : Neophyte + (1|Family/Genus), general)
  summary(MDL_LMA)
  MDL_freq<- lmer(impact ~ -1+Neophyte +Freq : Neophyte + (1|Family/Genus), general)
  summary(MDL_SM)
  MDL_freq_herb<- lmer(impact ~-1+Neophyte + Freq_herb : Neophyte + (1|Family/Genus), general)
  summary(MDL_SSD)
  MDL_soil<- lmer(impact ~-1+Neophyte + Soil : Neophyte + (1|Family/Genus), general)
  summary(MDL_N)
  MDL_mow<- lmer(impact ~-1+Neophyte + Mow : Neophyte + (1|Family/Genus), general)
  summary(MDL_N)
  MDL_gras<- lmer(impact ~-1+Neophyte + Gras : Neophyte + (1|Family/Genus), general)
  summary(MDL_N)
  
  
  
  # Plot
  visreg::visreg(MDL_sev,"Sev", by="Neophyte")
  visreg::visreg(MDL_sev_herb,"Sev_herb", by="Neophyte")
  visreg::visreg(MDL_freq,"Freq", by="Neophyte")
  visreg::visreg(MDL_freq_herb,"Freq_herb", by="Neophyte")
  visreg::visreg(MDL_soil,"Soil", by="Neophyte")
  visreg::visreg(MDL_mow,"Mow", by="Neophyte")
  visreg::visreg(MDL_gras,"Gras", by="Neophyte")
  ```
  
  
  
  ## 4.6 Effect all
  ```{r}
  # split analysis between woody and non-woody species?
  nw <- subset(general, woody != "woody",drop=T)
  w <- subset(general, woody == "woody", drop=T)
  alien <- subset(general, Neophyte != "native", drop=T)
  
  # scale predictors
  # Model for all traits together
  MDL<- lmer(impact ~-1+ SSD: Neophyte+ LMA: Neophyte+ H: Neophyte+ LA: Neophyte+ 
               SM: Neophyte+ N: Neophyte+ Neophyte+ woody +(1|Family/Genus), 
             general)
  summary(MDL)
  
  # Model for all traits together
  MDL<- lmer(impact ~-1+SM+SSD+LMA+H+LA+N+woody+Neophyte+(1|Family/Genus), 
             general)
  summary(MDL)
  
  impact[, c(20:25)] <- log(impact[, c(20:25)])
  impact[, c(20:25)] <- scale(impact[, c(20:25)])
  
  MDL<- lmer(RelDiff ~-1+SM+SSD+LMA+H+LA+N+woody+Neophyte+(1|Family/Genus), 
             impact[impact$class=="70%-100%",])
  summary(MDL)
  
  sum(!is.na(general$LA))
  sum(!is.na(general$LMA))
  sum(!is.na(general$N))
  sum(!is.na(general$SM))
  sum(!is.na(general$SSD))
  sum(!is.na(general$H))
  
  plot(DHARMa::simulateResiduals(MDL))
  
  # plot
  test_plot <- plot_model(MDL, type = "re", facet.grid=FALSE) 
  plot <- ggarrange(test_plot[[1]], test_plot[[2]], nrow=1, ncol=2, labels= c("a","b"),
                    font.label = list(size = 12))
  plot
  #ggsave("Random_effects.jpg", plot= plot, width = 10, height = 30)
  
  # check DHARMa
  simulation <- DHARMa::simulateResiduals(MDL)
  plot(simulation)
  ```
  
  
  
  ```{r}
  # Example interaction data
  # Define a sequence of SSD values over its range
  ssd_seq <- seq(min(general$SSD, na.rm = TRUE), max(general$SSD, na.rm = TRUE), length.out = 100)
  
  # Create a new dataset for predictions
  new_data <- expand.grid(
    SSD = ssd_seq,
    Neophyte = unique(general$Neophyte),
    LMA= mean(general$LMA, na.rm = TRUE),
    SM= mean(general$SM, na.rm = TRUE),
    H= mean(general$H, na.rm = TRUE),
    N= mean(general$N, na.rm = TRUE),
    LA= mean(general$LA, na.rm = TRUE), 
    Genus = DescTools::Mode(general$Genus, na.rm=T),
    Family =  unique((general$Family[general$Genus == DescTools::Mode(general$Genus, na.rm=T) & !is.na(general$Family)]))
  )
  
  # Add predictions to the dataset
  new_data$predicted <- predict(MDL, newdata = new_data, re.form = NA)  # Exclude random effects
  
  library(lme4)
  class(MDL)
  pred<- predict(MDL, newdata = new_data)
  ?predict
  new_data$lower <- pred[, "lwr"]
  new_data$upper <- pred[, "upr"]
  
  # Plot with confidence intervals
  ggplot(new_data, aes(x = SSD, y = predicted, color = Neophyte, group = Neophyte)) +
    geom_line(size = 1) + 
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = Neophyte), alpha = 0.2, color = NA) +
    labs(
      x = "SSD (Continuous)",
      y = "Predicted Impact",
      color = "Neophyte",
      fill = "Neophyte",
      title = "Interaction Effect: SSD x Neophyte with Confidence Intervals"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  ```
  
  
  ```{r}
  install.packages("marginaleffects")
  library(marginaleffects)
  newdata <- datagrid(newdata= general, by= (SSD = rep(seq(min(general$SSD, na.rm=T), 
                                                           max(general$SSD, na.rm=T), length=10), 3),
                                             Neophyte = c(rep("intra",10), rep("extra",10), rep("native",10))))
  colnames(general)
  predictions(MDL, newdata)
  predi
  
  
  library(effects)
  plot(effects::Effect(c("SSD","Neophyte"), MDL, se=T, confint=T, confidence.level= 0.05), multiline=F)
  ```
  
  
  ```{r}
  #
  library(ggeffects)
  plot(ggpredict(MDL,terms= c("LA", "Neophyte") ))
  
  get_model_data(MDL, type = "pred", terms= c("SSD", "Neophyte"), pred.type="re")
  
  plot_model(MDL, type = "pred", facet.grid=T, terms= c("SSD", "Neophyte"), pred.type="re", ci.lvl=0.95, show.data=T, alpha=0.1) 
  plot_model(MDL, type = "pred", facet.grid=T, terms= c("LMA", "Neophyte"), pred.type="re", ci_level=NA, show.data=T) 
  plot_model(MDL, type = "pred", facet.grid=T, terms= c("LA", "Neophyte"), pred.type="re", ci_level=NA, show.data=T) 
  plot_model(MDL, type = "pred", facet.grid=T, terms= c("N", "Neophyte"), pred.type="re", ci_level=NA, show.data=T) 
  plot_model(MDL, type = "pred", facet.grid=T, terms= c("H", "Neophyte"), pred.type="re", ci_level=NA, show.data=T) 
  plot_model(MDL, type = "pred", facet.grid=T, terms= c("SM", "Neophyte"), pred.type="re", ci_level=NA, show.data=T) 
  
  visreg::visreg(MDL, "Neophyte")
  ```
  
  
  ## 4.7 PCA
  ```{r}
  # look at number of traits available per species
  eva_names$traits <- Diaz$`Number of traits with values`[match(eva_names$name, Diaz$name)]
  eva_names$traits[is.na(eva_names$traits)] <- 0
  
  # look at available data
  hist(eva_names$traits[eva_names$name %in% general$taxa])
  
  # generate dataset all data available
  complete <- general[general$number=="6" & !is.na(general$number),]
  
  complete_nw <- subset(complete, woody != "woody", drop=T)
  complete_w <- subset(complete, woody== "woody", drop=T)
  
  MDL<- lmer(impact ~-1+ SSD: Neophyte+ LMA: Neophyte+ H: Neophyte+ LA: Neophyte+ SM: Neophyte+ 
               N: Neophyte+ Neophyte +(1|Family/Genus), complete_nw)
  summary(MDL)
  
  # test differences PCA alien native
  complete_nw_intra <- subset(complete_nw, Neophyte == "intra", drop=T)
  complete_nw_native <- subset(complete_nw, Neophyte == "native", drop=T)
  complete_nw_extra <- subset(complete_nw, Neophyte == "extra", drop=T)
  
  # test PCA
  PCA <- prcomp(complete_nw[,c(3:8)], scale=F)
  eig.val <- get_eigenvalue(PCA)
  eig.val
  biplot(PCA)
  PCA
  
  # plot PCA
  library(ggfortify)
  autoplot(PCA, data = complete_nw,
           loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3)+
    geom_point(aes(size=Neophyte, colour= Neophyte))
  
  # test location alien species 
  test <- cbind(complete_nw[,2], as.data.frame(PCA$x[, 1:2]))
  test$Neophyte <- as.factor(test$Neophyte)
  result<-welch_anova_test(test, PC1 ~ Neophyte)
  result
  games_howell_test(data= test, PC1 ~ Neophyte)  
  ```
  
  
  # 5 Comparison
  ## 5.1 Status
  ```{r}
  # Data on which species are neophytes
  native_intra_analysis=F
  if(native_intra_analysis){
    species_country_status<- read_csv("country_species_ESy.csv", show_col_types = FALSE)
  } else{
    species_country_status<- read_csv("country_species_ESy.csv", show_col_types = FALSE)
    species_country_status$Neophyte[species_country_status$Neophyte=="native_intra"] <- "native"
  }
  
  # make list unique species and bind with traits
  species <- unique(species_country_status[,c(5,9)])
  species <- cbind(species, eva[match(species$name, eva$name), c(33:40)])
  
  # look at mean values
  traits <- species |> group_by(Neophyte) |> 
    summarise(n=n(), 
              LMA = mean(`LMA (g/m2)`, na.rm=T),
              SSD = mean(`SSD combined (mg/mm3)`, na.rm=T),
              LA= mean(`Leaf area (mm2)`, na.rm=T),
              N = mean(`Nmass (mg/g)`, na.rm=T),
              SM= mean(`Diaspore mass (mg)`, na.rm=T),
              H= mean(`Plant height (m)`, na.rm=T),
              rel = sum(`Growth Form`=="tree", na.rm=T)/n)
  
  # scale data
  species[, c(5:10)] <- scale(log(species[, c(5:10)]))
  
  # change levels
  species <- species %>% mutate(Neophyte = factor(Neophyte, 
                                                  levels = c("native", "intra","extra"),
                                                  labels = c("native in the country", 
                                                             "intra-European neophyte", 
                                                             "extra-European neophyte")))
  
  not_tree <- species[!species$`Growth Form`=="tree" & !is.na(species$`Growth Form`),]
  
  traits_not_tree <- not_tree |> group_by(Neophyte) |> 
    summarise(n=n(), 
              LMA = mean(`LMA (g/m2)`, na.rm=T),
              SSD = mean(`SSD combined (mg/mm3)`, na.rm=T),
              LA= mean(`Leaf area (mm2)`, na.rm=T),
              N = mean(`Nmass (mg/g)`, na.rm=T),
              SM= mean(`Diaspore mass (mg)`, na.rm=T),
              H= mean(`Plant height (m)`, na.rm=T),
              rel = sum(`Growth Form`=="tree", na.rm=T)/n)
  
  ggplot(not_tree,
         aes(x= Neophyte, y=`LMA (g/m2)`, fill=Neophyte, 
             color= Neophyte))+
    geom_violin(alpha=0.5, scale="width")+
    geom_boxplot(width= 0.25, alpha=0.8, fill="white")+
    theme_pubr()+
    stat_summary(fun= "mean",
                 geom = "point", aes(group= Neophyte), size=3)+
    scale_colour_manual(values=c("#1E88E5", "#FFC107", "#004D40"))+
    scale_fill_manual(values = c("#1E88E5", "#FFC107", "#004D40")) +
    theme(legend.position = "none")
  
  # test differences and significance
  result<-welch_anova_test(not_tree, `LMA (g/m2)` ~ Neophyte )
  result
  games_howell_test(data= not_tree, `LMA (g/m2)` ~ Neophyte)
  
  result<-welch_anova_test(not_tree, `Nmass (mg/g)` ~ Neophyte )
  result
  
  result<-welch_anova_test(species, `Diaspore mass (mg)` ~ Neophyte )
  result
  games_howell_test(data= species, `Diaspore mass (mg)`~ Neophyte)
  
  result<-welch_anova_test(species, `Leaf area (mm2)` ~ Neophyte )
  result
  games_howell_test(data= species, `Leaf area (mm2)`~ Neophyte)
  
  
  result<-welch_anova_test(species[!species$`Growth Form`=="tree",], `Plant height (m)` ~ Neophyte )
  result
  games_howell_test(data= species[!species$`Growth Form`=="tree",], `Plant height (m)` ~ Neophyte)
  ```
  
  
  
  # 6 TREE
  ```{r}
  library(V.PhyloMaker2)
  phylo <- data.frame(name= general$taxa, genus= general$impact, family= general$Family)
  tree <- phylo.maker(sp.list= phylo, tree= GBOTB.extended.TPL, nodes= nodes.info.1.TPL, scenarios = "S3")
  #write.tree(tree$scenario.3, "sample.tre")
  ```
  
  
  Plot
  ```{r}
  library(phylogram)
  library(ggtree)
  
  ggtree(tree2$scenario.1, layout="circular", ladderize = FALSE)+
    geom_nodepoint() + geom_tiplab(hjust = -.1)
  ```
  
  ###### GAM PLOTS ####
  
  ```{r}
  # Load required libraries
  library(mgcv)
  library(gratia)
  library(ggplot2)
  
  # Choose the term you want to plot
  term_to_plot <- "s(EIVEresT)"  # change as needed
  
  draw(MDL, residuals = TRUE, select= "s(EIVEresT)")
  
  # Create partial residual plot on the RESPONSE scale
  partial_plot <- draw(
    mod,
    select = term_to_plot,
    residuals = TRUE,
    transform = exp  # for negative binomial to response scale
  )
  
  # Show the plot
  print(partial_plot)
  
  sm <- smooth_estimates(MDL) |>
    add_confint()
  sm
  
  # add partial residuals to data
  eg1 <- fullPlotData |>
    add_partial_residuals(MDL)
  names(eg1)
  eg1 <- as.data.frame(eg1)
  dim(eg1)
  p_sx2 <- sm |> filter(.smooth == paste("s(", term_to_plot,")",sep="") |> ggplot() + 
                          geom_rug(aes(x = term_to_plot), data = eg1, length = grid::unit(0.02, "npc"))+
                          geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci, x = x2), alpha = 0.2))
  
  geom_point(aes(x = x2, y = `s(x2)`),
             data = eg1, cex = 1.5, colour = "steelblue3"
  ) +
    geom_line(aes(x = x2, y = .estimate), lwd = 1.2) +
    labs(y = "Partial effect", title = exp())
  p_sx2
  
  plot(MDL, trans=  exp)
  
  max(MDL$family$linkinv(preds$fit))
  
  
  names(coef(MDL))
  install.packages("gratia")
  library(gratia)
  
  # Example for a single variable `var_name`
  var_name <- "EIVEresT"  # change this in a loop
  label <- "EIVE T"
  col <- "#1E88E5"
  
  # Get partial residuals
  partial_resids <- partial_residuals(MDL, select = var_name, partial_match=T)
  
  # Create grid for predictions
  newdata <- with(fullPlotData, data.frame(!!var_name <- seq(min(get(var_name)), max(get(var_name)), length.out = 100)))
  # Fill in other required covariates with their means or representative values
  other_vars <- setdiff(names(coef(MDL)), var_name)
  for (v in other_vars) {
    if (!v %in% names(newdata)) newdata[[v]] <- mean(all[[v]], na.rm = TRUE)
  }
  
  
  
  # Get model predictions (on link scale)
  newdata$fit <- predict(MDL, newdata = newdata, type = "link", se.fit = TRUE)$fit
  newdata$se <- predict(MDL, newdata = newdata, type = "link", se.fit = TRUE)$se.fit
  newdata$lower <- newdata$fit - 1.96 * newdata$se
  newdata$upper <- newdata$fit + 1.96 * newdata$se
  
  # Plot
  ggplot() +
    geom_point(data = partial_resids, aes_string(x = var_name, y = ".partial.resid"), alpha = 0.4) +
    geom_line(data = newdata, aes_string(x = var_name, y = "fit"), color = col, linewidth = 1.2) +
    geom_ribbon(data = newdata, aes_string(x = var_name, ymin = "lower", ymax = "upper"), alpha = 0.2, fill = col) +
    labs(x = label, y = "Partial effect (link scale)") +
    theme_pubr()
  
  new_data <- as.data.frame(matrix(NA,nrow=1000, ncol= 17))
  colnames(new_data) <- c("Area","EIVEresT","EIVEresN","EIVEresM","EIVEresR","EIVEresL","DistSeverity.sqrt",        
                          "Soil.Disturbance.sqrt","Grazing.Pressure.sqrt","Mowing.Frequency.sqrt","hfp","elev","chelsaP","DI_extra", 
                          "DI_intra", "Latitude","Longitude")
  new_data[, colnames(new_data)==var[i,1]] <- seq(min(fullPlotData[, colnames(fullPlotData)== var[i,1]]), 
                                                  max(fullPlotData[, colnames(fullPlotData)== var[i,1]]), length.out = 1000)
  others <- strsplit(as.character(MDL$pred.formula)[2],"[/+]")[[1]]
  others <- gsub(" ","",others)
  others <- others[!others %in% var[i,1]]
  
  for(v in 1:16){
    index <- others[v]
    new_data[, colnames(new_data)==index] <- rep(mean(fullPlotData[[index]], na.rm=T),times = 1000)
  }
  
  
  preds <- predict(MDL, type="terms", newdata= new_data, se.fit=T)
  
  
  plot_data <- as.data.frame(matrix(NA, nrow=1000, ncol=1))
  plot_data[,1] <- new_data[, colnames(new_data)==var[i,1]]
  plot_data$pred <- preds$fit[, colnames(preds$fit)== paste("s(",var[i,1], ")", sep="")]
  plot_data$up <- plot_data$pred +1.96*preds$se.fit[, colnames(preds$se.fit)== paste("s(",var[i,1], ")", sep="")]
  plot_data$down <- plot_data$pred -1.96*preds$se.fit[, colnames(preds$se.fit)== paste("s(",var[i,1], ")", sep="")]   
  
  plot_end <- ggplot(plot_data, aes(V1, y= pred, color=col, fill= col))+
    geom_line()+
    geom_ribbon(aes(ymin= down, ymax=up), alpha=0.2)
  plot_end
  
  plot(plot_data$V1, plot_data$pred)
  
  plot(plot_pred, show_residuals=T)
  
  test <- residualize_over_grid(plot_pred, MDL, protect_names = T)
  
  
  ```
  
  