
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