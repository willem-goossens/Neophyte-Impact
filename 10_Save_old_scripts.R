
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
  
  