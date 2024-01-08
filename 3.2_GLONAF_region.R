rm(list=ls())

##### LOAD ####
library(readr)
library(ggplot2)
library(sf)
#library(rgeos)
library(sp)
library(units)
library(tibble)
library(tidyverse)
library(sp)
library(units)
library(lwgeom)
library(terra)
library(diffdf)
library(doParallel)
library(arsenal)

##### DATA ####
###### header ######
header <- read_delim("../EVA Data/171_NeophyteInvasions20230216_notJUICE_header.csv", "\t")
fast <- T
if(fast) {
  header <- header[runif(length(header$PlotObservationID)) > 0.9999,]
}

header <- header[!(is.na(header$Latitude) | is.na(header$Longitude)),]

plotLocations <- st_as_sf(header, coords = c("Longitude","Latitude"), remove = FALSE)
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")

###### glonaf ######
glonafRegions <- read_sf("../GloNAF_Shapefile", "regions2")
glonafRegions <- st_transform(glonafRegions, CRS("+proj=longlat +datum=WGS84"))

glonafRegionList<- read.csv("../GloNAF_Shapefile/Region_GloNAF_vanKleunenetal2018Ecology.csv")
glonafRegionList<- glonafRegionList[glonafRegionList$tdwg1_name=="Europe", ]  
unique(glonafRegionList$tdwg4_name)
newGlonaf<- subset(glonafRegions, glonafRegions$OBJIDsic %in% glonafRegionList$OBJIDsic)
glonafRegions<- newGlonaf
all.equal(sort(glonafRegionList$OBJIDsic), sort(glonafRegions$OBJIDsic))
newGlonaf<- left_join(newGlonaf, glonafRegionList, by= "OBJIDsic")

#(st_geometry(glonafRegions), reset=F)
#plot(st_geometry(glonafRegions[24,]), col="red", add=T)
#plot(st_geometry(plotLocations), add=T)

###### joined ######
sf_use_s2(T)
joinedData<-st_join( plotLocations, newGlonaf[c(1:23,25:85),],join=st_within)
# some data is duplicated because it is located on the edge of two countries and glonaf shapefile is not very accurate
setdiff(joinedData$PlotObservationID, plotLocations$PlotObservationID)
any(is.na(joinedData$PlotObservationID))
any(is.na(plotLocations$PlotObservationID))
joinedData$PlotObservationID[duplicated(joinedData$PlotObservationID)]
# Here we just assign it but has to be improved
joinedData<-joinedData[!duplicated(joinedData$PlotObservationID),]

# Northern Ireland cannot be merged using s2
sf_use_s2(F)
res <- st_join(plotLocations, newGlonaf[24,])
any(!is.na(res$OBJIDsic))
index<-which(!is.na(res$OBJIDsic))
# We cannot merge Northern Ireland
joinedData[index,] <- res[index,]
Northern_Ireland<- res[!is.na(res$OBJIDsic),]

#pi= st_intersection(glonafRegions[glonafRegions$OBJIDsic=="457",],glonafRegions[glonafRegions$OBJIDsic=="1558",])
#plot(st_geometry(glonafRegions[glonafRegions$OBJIDsic=="457",]), col='blue')
#plot(st_geometry(glonafRegions[glonafRegions$OBJIDsic=="1558",]), add=T, col="yellow")
#plot(pi$geometry, add=T, col="red")

percentNotAssignedPlots <- sum(is.na(joinedData$code))/length(joinedData$code)*100
percentNotAssignedPlots
length(joinedData$code)

all.equal(joinedData$PlotObservationID, plotLocations$PlotObservationID)


##### NOT ASSIGNED #####
###### original #####
#distanceThreshold <- set_units(10000,m)

#remainingPlots <- plotLocations[is.na(joinedData$code),]
#remainingPlots<- remainingPlots[, c("PlotObservationID","Country", "geometry")]

# Compute bounding boxes for the regions, this allows to decide faster, if a plot is outside the distanceThreshold from a region
#boundingBoxes <- glonafRegions

#for(i in 1:length(boundingBoxes$OBJIDsic)) {
#  boundingBoxes$geometry[i] <- st_as_sfc(st_bbox(boundingBoxes$geometry[i]))
#}

#remainingPlots$code <- NA
#remainingPlots$Distance <- -1

#begin<-Sys.time()  
#for(i in 1:length(remainingPlots$PlotObservationID)) {
#  region <- NA
#  distance <- set_units(1000000,m)
#
#for(j in 1:length(glonafRegions$OBJIDsic)) {
#    distBoundingBox <- st_distance(boundingBoxes[j,], remainingPlots[i,]) #Calculate the distance
#   # Only if the distance to the bounding box is smalle than distanceThreshold we need to compute the real distance to the region
#   if(distBoundingBox < distanceThreshold){
#     dist <- st_distance(glonafRegions[j,], remainingPlots[i,])
#     if(dist < distanceThreshold & dist < distance) {
#             region <- glonafRegions$OBJIDsic[j]
#       distance <- dist
#     }
#   }
# }
# Print the progress
#  message(i/length(remainingPlots$PlotObservationID)*100, "% ", remainingPlots$Country[i], " ", region)
#  remainingPlots$code[i] <- region
# remainingPlots$Distance[i] <- distance
#}
#end<-Sys.time()
#round(end-begin, 2)

#remainingPlots2<- remainingPlots

#allPlotsWithRegion <- joinedData[,c("PlotObservationID", "code", "Country")]
#allPlotsWithRegion[is.na(allPlotsWithRegion$code),] <- remainingPlots[,c("PlotObservationID", "code", "Country")]

#x<-allPlotsWithRegion[is.na(allPlotsWithRegion$code),]

###### parallel #####
# We will run it parallel
parallel::detectCores()
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()
foreach::getDoParWorkers()

distanceThreshold <- set_units(10000,m)

remainingPlots <- plotLocations[is.na(joinedData$code),]

# make smaller to make calculation quicker
remainingPlots<- remainingPlots[, c("PlotObservationID","Country", "geometry")]

# Compute bounding boxes for the regions, this allows to decide faster, if a plot is outside the distanceThreshold from a region
boundingBoxes <- glonafRegions

for(i in 1:length(boundingBoxes$OBJIDsic)) {
  boundingBoxes$geometry[i] <- st_as_sfc(st_bbox(boundingBoxes$geometry[i]))
}

remainingPlots$code <- NA
remainingPlots$Distance <- -1

length(remainingPlots$PlotObservationID)
begin<-Sys.time() 
remainingPlots[,4:5]<-foreach(i = 1:length(remainingPlots$PlotObservationID), .combine='rbind', .packages=c("dplyr","mgcv", "sf","units")) %dopar% {
  region <- NA
  distance <- set_units(1000000,m)
  
  for(j in 1:length(boundingBoxes$OBJIDsic)) {
    sf_use_s2(F)
    distBoundingBox <- st_distance(boundingBoxes[j,], remainingPlots[i,]) #Calculate the distance
    # Only if the distance to the bounding box is smalle than distanceThreshold we need to compute the real distance to the region
    if(distBoundingBox < distanceThreshold){
      dist <- st_distance(glonafRegions[j,], remainingPlots[i,])
      if(dist < distanceThreshold & dist < distance) {
        region <- glonafRegions$OBJIDsic[j]
        distance <- dist
      }
    }
  }
  # Print the progress
  message(i/length(remainingPlots$PlotObservationID)*100, "% ", remainingPlots$Country[i], " ", region)
  remainingPlots$code[i] <- c(region, distance)
}
end<-Sys.time()
round(end-begin, 2)  

remainingPlots$tdwg4_name<- NA
for (i in 1: length(unique(glonafRegionList$OBJIDsic))) {
  index<- unique(glonafRegionList$OBJIDsic)[i]
  remainingPlots$tdwg4_name[remainingPlots$code==index]<- glonafRegionList$tdwg4_name[glonafRegionList$OBJIDsic== index]
}

remainingPlots
#print(setdiff(remainingPlots$code, remainingPlots2$code))

parallel::stopCluster(cl = my.cluster)

# merge the plots which lie within a region with those which are in proximity and write the result to file
allPlotsWithRegion <- joinedData[,c("PlotObservationID", "tdwg4_name", "Country")]
remainingPlots$code<- as.character(remainingPlots$code)
allPlotsWithRegion[is.na(allPlotsWithRegion$tdwg4_name),] <- remainingPlots[,c("PlotObservationID", "tdwg4_name", "Country")]
any(is.na(allPlotsWithRegion$tdwg4_name))

all.equal(sum(is.na(allPlotsWithRegion$tdwg4_name)), sum(is.na(remainingPlots$tdwg4_name)))

colnames(allPlotsWithRegion)<- c("PlotObservationID", "Region", "Country", "geometry")

#comparedf(st_drop_geometry(remainingPlots), st_drop_geometry(remainingPlots2))

###### assign #####
#write.csv(allPlotsWithRegion,"allPlotsWithRegion_small.csv")

compareCountryLabelWithRegion <- function(country, region) {
  message("Number of plots in region: ", sum(allPlotsWithRegion$Region == region, na.rm = T))
  message("Number of plots in country: ", sum(allPlotsWithRegion$Country == country, na.rm = T))
}

# Compare the number of plots in a region with the number of plots which are in the "same" country, does it make sense?
compareCountryLabelWithRegion("Albania", "Albania")
compareCountryLabelWithRegion("Austria", "Austria")
compareCountryLabelWithRegion("Netherlands", "Netherlands")
compareCountryLabelWithRegion("Italy", "Italy")
compareCountryLabelWithRegion("Sicily", "Sicily")
compareCountryLabelWithRegion("Italy", "Sardinia")

unique(allPlotsWithRegion$Region)
unique(allPlotsWithRegion$Country)
