##### 1 START #####
rm(list=ls())
# load packages
library(terra)
library(tidyterra)
library(sf)
library(geodata)
library(sp)
library(readr)
library(ggplot2)
library(doParallel)
library(foreach)
library(parallel)
library(graphics)
library(lubridate)
library(climenv)
library(units)

##### 2 LOAD ####
###### 2.1 HFP #####
# load Human Foot Print data and check crs
hfp2009<- rast("../Disturbance Data/HFP2009.tif")
hfp1993<- rast("../Disturbance Data/HFP1993.tif")
cat(crs(hfp2009))

###### 2.2 MED regions #####
# load medregions data and check crs
medRegions <- read_sf("../Europe-regions-shapefiles-2023", "Emed_regions")
medRegions <- st_transform(medRegions, CRS("+proj=longlat +datum=WGS84"))
cat(crs(medRegions))

###### 2.3 Plots #####
# load full plot data
fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
fullPlotData<- fullPlotData[,c(1:5)]

# reduce size
fast <- F
if(fast) {
  fullPlotData <- fullPlotData[runif(length(fullPlotData$PlotObservationID)) > 0.8,]
}

# Get plot locations and give geometry
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")


##### 3 ANALYSIS ####
###### 3.1 Prepare #####
# crop human foot print data using a bounding box around europe
mollweide_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"
europe_bbox <- st_bbox(c(xmin = -3000000, xmax = 6000000, ymin = 3500000, ymax = 10000000), crs = mollweide_crs)
europe_polygon <- st_as_sfc(europe_bbox)
hfp2009<- crop(hfp2009, europe_polygon)
hfp1993<- crop(hfp1993, europe_polygon)

# change projection to that of medRegion
hfp2009<- terra::project(hfp2009, "+proj=longlat +datum=WGS84")
cat(crs(hfp2009))
hfp1993<- terra::project(hfp1993, "+proj=longlat +datum=WGS84")

# crop medregions from human footprint map
hfp2009<- mask(hfp2009, medRegions)
hfp2009<- crop(hfp2009, medRegions)
hfp1993<- mask(hfp1993, medRegions)
hfp1993<- crop(hfp1993, medRegions)

# plot data to check
plot(hfp2009)   # QUESTION: why is part lower masked?
plot(hfp1993) 


###### 3.2 HFP 2009 #####
# Extract the values for HFP 2009
test<- terra::extract(hfp2009, plotLocations[,c(2:4)])
test<- test[,2]
# Combine 
plotLocations<- cbind(plotLocations, test)

# Show
plotting=F
if(plotting){
ggplot() +   
  geom_spatraster(data = hfp2009) + 
  geom_sf(data = plotLocations, color = "black", size = 1) + 
  coord_sf()
}

# We are not able to extract the value for all sites so we assign the nearest value (TO BE DISCUSSED)  


x<- read_csv("remainingPlotsHFP2009.csv")
plotLocations$test[is.na(plotLocations$test)] <- x$test[match(plotLocations$PlotObservationID[is.na(plotLocations$test)], x$PlotObservationID)]

remainingPlots<- plotLocations[is.na(plotLocations$test),]
# I was not able to do this with foreach so takes about 1 hour to calculate
for(i in 1:nrow(remainingPlots)){
  # Buffer point locations by desired distance
  ptBuff <- st_buffer(remainingPlots[i,], dist = 50000)
  # Crop raster to buffered point
  inbuff <- crop(hfp2009, vect(ptBuff))
  # Convert to points
  dat <- as.points(inbuff)
  dat <- st_as_sf(dat)
  # If all values are NA the dataframe will be empty, move to next
  if(nrow(dat) == 0){
    next
  }
  ptdist <- st_distance(remainingPlots[i,], dat)
  dat$pdistance <- as.numeric(ptdist)
  m <- as.numeric(min(ptdist))
  val <- dat %>% 
    dplyr::filter(pdistance == m) %>% 
    select(HFP2009) %>% 
    st_drop_geometry()
  val <- as.numeric(val)
  remainingPlots[i,]$test <- val
}

# assign the new values to the previous NAs  
#plotLocations[is.na(plotLocations$test),] <- remainingPlots
#write_csv(plotLocations, "remainingPlotsHFP2009.csv")

# Make dataset with all plots that are still not assigned --> probably due to the raster nature --> maybe make the distance larger
remainingPlots<- plotLocations[is.na(plotLocations$hfp2009),]
# Plot these remaining not extracted plots --> all close to the coast
if(plotting){
ggplot() +   
  geom_spatraster(data = hfp2009) + 
  geom_sf(data = remainingPlots, color = "black", size = 1) + 
  coord_sf()+theme_minimal()
}

# prepare for saving
saving= T
if(saving){
  st_geometry(plotLocations) <- NULL
  colnames(plotLocations)[6]<- "hfp2009"
  fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
  fullPlotData<- left_join(fullPlotData, plotLocations)
  write_csv(fullPlotData, "fullPlotData_ESy.csv")
}

###### 3.3 HFP 1993 #####
# Get plot locations and give geometry
if(fast) {
  fullPlotData <- fullPlotData[runif(length(fullPlotData$PlotObservationID)) > 0.99,]
}


plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5)]
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")

# Extract the values for HFP 1993
test<- terra::extract(hfp1993, plotLocations[,c(2:4)])
test<- test[,2]
# Combine 
plotLocations<- cbind(plotLocations, test)

# Show
plotting=F
if(plotting){
  ggplot() +   
    geom_spatraster(data = hfp1993) + 
    geom_sf(data = plotLocations, color = "black", size = 1) + 
    coord_sf()
}

x <- read_csv("remainingPlotsHFP1993.csv")
plotLocations$test[is.na(plotLocations$test)] <- x$test[match(plotLocations$PlotObservationID[is.na(plotLocations$test)], x$PlotObservationID)]

# We are not able to extract the value for all sites so we assign the nearest value (TO BE DISCUSSED)  
remainingPlots<- plotLocations[is.na(plotLocations$test),]


# I was not able to do this with foreach so takes about 1 hour to calculate
for(i in 1:nrow(remainingPlots)){
  # Buffer point locations by desired distance
  ptBuff <- st_buffer(remainingPlots[i,], dist = 50000)
  # Crop raster to buffered point
  inbuff <- crop(hfp1993, vect(ptBuff))
  # Convert to points
  dat <- as.points(inbuff)
  dat <- st_as_sf(dat)
  # If all values are NA the dataframe will be empty, move to next
  if(nrow(dat) == 0){
    next
  }
  ptdist <- st_distance(remainingPlots[i,], dat)
  dat$pdistance <- as.numeric(ptdist)
  m <- as.numeric(min(ptdist))
  val <- dat %>% 
    dplyr::filter(pdistance == m) %>% 
    select(HFP1993) %>% 
    st_drop_geometry()
  val <- as.numeric(val)
  remainingPlots[i,]$test <- val
}

# assign the new values to the previous NAs  
plotLocations[is.na(plotLocations$test),] <- remainingPlots

write_csv(plotLocations, "remainingPlotsHFP1993.csv")

# Make dataset with all plots that are still not assigned --> probably due to the raster nature --> maybe make the distance larger
remainingPlots<- plotLocations[is.na(plotLocations$test),]
# Plot these remaining not extracted plots --> all close to the coast
if(plotting){
  ggplot() +   
    geom_spatraster(data = hfp1993) + 
    geom_sf(data = remainingPlots, color = "black", size = 1) + 
    coord_sf()+theme_minimal()
}


# prepare for saving
saving= T
if(saving){
  st_geometry(plotLocations) <- NULL
  fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
  colnames(plotLocations)[6]<- "hfp1993"
  fullPlotData<- left_join(fullPlotData, plotLocations)
  write_csv(fullPlotData, "fullPlotData_ESy.csv")
  #write_csv(fullPlotData, "fullPlotData_euro.csv")
}
# Note that the lengths are not of the same length --> strange?
# for safety reasons also did it in _euro

#fullPlotData <- read_csv("fullPlotData_ESy.csv")

###### 3.4 Time #####
# Load dataset again if not wanting to do all previous work
# fullPlotData<- read_csv("fullPlotData_euro.csv", show_col_types = FALSE)

# Check if there are any NAs --> not here but yes later --> are dates 00.00.0000 --> remove later!
anyNA(fullPlotData$Date)
# Make dates
fullPlotData$Date<- as.Date(fullPlotData$Date, format = "%d.%m.%Y")

# Assign extracted data to correct plot
fullPlotData$hfp<- ifelse(fullPlotData$Date > dmy("01012001"), fullPlotData$hfp2009, fullPlotData$hfp1993)
nrow(fullPlotData[(is.na(fullPlotData$hfp)),])

saving= T
if(saving){
  write_csv(fullPlotData, "fullPlotData_ESy.csv")
}

# TO DO: think about too old plots (is 1870s still represented by the level of human pressure 100 years later)


##### 4 CLIMATE ####
###### 4.1 Prepare ####
# Load data
chelsaT<- rast("../EIVE Data/bio1.tif")
chelsaP<- rast("../EIVE Data/bio12.tif")
elev<- rast("../EIVE Data/Elev.tif")

# create box around Europe
europe_bbox <- st_bbox(c(xmin = -30, xmax = 80, ymin = 25, ymax = 90))
europe_polygon <- st_as_sfc(europe_bbox)

# crop this box and reproject
chelsaT<- crop(chelsaT, europe_polygon)
chelsaT<- terra::project(chelsaT, "+proj=longlat +datum=WGS84")
chelsaP<- crop(chelsaP, europe_polygon)
chelsaP<- terra::project(chelsaP, "+proj=longlat +datum=WGS84")
elev<- crop(elev, europe_polygon)
elev<- terra::project(elev, "+proj=longlat +datum=WGS84")
# do we mask or not? there is something to say for not doing it (as we then really use this coordinate present in the plotdata file) but we do it for
# the other data, so also using a buffer might be a better idea
#chelsaT<- terra::mask(chelsaT, medRegions) 
#chelsaP<- terra::mask(chelsaP, medRegions) 
#elev<- terra::mask(elev, medRegions) 

# crop medRegions and plot
chelsaT<- terra::crop(chelsaT, medRegions)
plot(chelsaT)
chelsaP<- terra::crop(chelsaP, medRegions)
plot(chelsaP)
elev<- terra::crop(elev, medRegions)
plot(elev)


###### 4.2 T mean #####
# Get plot locations and give geometry
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5)]
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")

# Extract the values for chelsa T
test<- terra::extract(chelsaT, plotLocations[,c(2:4)])
test<- test[,2]
# Combine 
plotLocations<- cbind(plotLocations, test)

# Show
plotting=F
if(plotting){
  ggplot() +   
    geom_spatraster(data = chelsaT) + 
    geom_sf(data = plotLocations, color = "black", size = 0.0001) + 
    coord_sf()+scale_fill_continuous(low="darkblue", high="darkred", 
                                     guide="colorbar",na.value="white")+
    theme_classic()
}

# We are not able to extract the value for all sites so we assign the nearest value (TO BE DISCUSSED)  
remainingPlots<- plotLocations[is.na(plotLocations$test),]


# I was not able to do this with foreach
for(i in 1:nrow(remainingPlots)){
  # Buffer point locations by desired distance
  ptBuff <- st_buffer(remainingPlots[i,], dist = 10000)
  # Crop raster to buffered point
  inbuff <- crop(chelsaT, vect(ptBuff))
  # Convert to points
  dat <- as.points(inbuff)
  dat <- st_as_sf(dat)
  # If all values are NA the dataframe will be empty, move to next
  if(nrow(dat) == 0){
    next
  }
  ptdist <- st_distance(remainingPlots[i,], dat)
  dat$pdistance <- as.numeric(ptdist)
  m <- as.numeric(min(ptdist))
  val <- dat %>% 
    dplyr::filter(pdistance == m) %>% 
    select(bio1) %>% 
    st_drop_geometry()
  val <- as.numeric(val)
  remainingPlots[i,]$test <- val
}

# assign the new values to the previous NAs  
plotLocations[is.na(plotLocations$test),] <- remainingPlots
#write_csv(plotLocations, "remainingPlots_T.csv")


# Make dataset with all plots that are still not assigned --> probably due to the raster nature --> maybe make the distance larger
remainingPlots<- plotLocations[is.na(plotLocations$test),]
# Plot these remaining not extracted plots --> all close to the coast
if(plotting){
  ggplot() +   
    geom_spatraster(data = chelsaT) + 
    geom_sf(data = remainingPlots, color = "black", size = 1) + 
    coord_sf()+theme_minimal()
}


# prepare for saving
saving= T
if(saving){
  st_geometry(plotLocations) <- NULL
  fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
  colnames(plotLocations)[6]<- "chelsaT"
  fullPlotData<- left_join(fullPlotData, plotLocations)
  write_csv(fullPlotData, "fullPlotData_ESy.csv")
}



###### 4.3 P mean #####
# Get plot locations and give geometry
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5)]
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")

# Extract the values for chelsa P
plotLocations <- st_as_sf(plotLocations, coords = c("Longitude","Latitude"), remove = FALSE)
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")

test<- terra::extract(chelsaP, plotLocations[,c(2:4)])
test<- test[,2]
# Combine 
plotLocations<- cbind(plotLocations, test)

# Show
plotting=F
if(plotting){
  ggplot() +   
    geom_spatraster(data = chelsaP, aes(fill=bio12)) + 
    geom_sf(data = plotLocations, color = "black", size = 0.001) + 
    coord_sf()+scale_fill_continuous(low="lightblue", high="darkblue", 
                                      guide="colorbar",na.value="white")+
    theme_classic()
}

# We are not able to extract the value for all sites so we assign the nearest value (TO BE DISCUSSED)  
remainingPlots<- plotLocations[is.na(plotLocations$test),]


# I was not able to do this with foreach so takes about 1 hour to calculate
for(i in 1:nrow(remainingPlots)){
  # Buffer point locations by desired distance
  ptBuff <- st_buffer(remainingPlots[i,], dist = 10000)
  # Crop raster to buffered point
  inbuff <- crop(chelsaP, vect(ptBuff))
  # Convert to points
  dat <- as.points(inbuff)
  dat <- st_as_sf(dat)
  # If all values are NA the dataframe will be empty, move to next
  if(nrow(dat) == 0){
    next
  }
  ptdist <- st_distance(remainingPlots[i,], dat)
  dat$pdistance <- as.numeric(ptdist)
  m <- as.numeric(min(ptdist))
  val <- dat %>% 
    dplyr::filter(pdistance == m) %>% 
    select(bio12) %>% 
    st_drop_geometry()
  val <- as.numeric(val)
  remainingPlots[i,]$test <- val
}

# assign the new values to the previous NAs  
plotLocations[is.na(plotLocations$test),] <- remainingPlots
#write_csv(plotLocations, "remainingPlots_P.csv")

# Make dataset with all plots that are still not assigned --> probably due to the raster nature --> maybe make the distance larger
remainingPlots<- plotLocations[is.na(plotLocations$test),]
# Plot these remaining not extracted plots --> all close to the coast
if(plotting){
  ggplot() +   
    geom_spatraster(data = chelsaP) + 
    geom_sf(data = remainingPlots, color = "black", size = 1) + 
    coord_sf()+theme_minimal()
}

# prepare for saving
saving= T
if(saving){
  st_geometry(plotLocations) <- NULL
  colnames(plotLocations)[6]<- "chelsaP"
  fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
  fullPlotData<- left_join(fullPlotData, plotLocations)
  write_csv(fullPlotData, "fullPlotData_ESy.csv")
}


###### 4.3 Elev mean #####
# Get plot locations and give geometry
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5)]
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")

# Extract the values for chelsa P
plotLocations <- st_as_sf(plotLocations, coords = c("Longitude","Latitude"), remove = FALSE)
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")

test<- terra::extract(elev, plotLocations[,c(2:4)])
test<- test[,2]
# Combine 
plotLocations<- cbind(plotLocations, test)

# Show
plotting=T
if(plotting){
  ggplot() +   
    geom_spatraster(data = elev, aes(fill=Elev)) + 
    #geom_sf(data = plotLocations, color = "black", size = 0.001) + 
    coord_sf()+scale_fill_continuous(low="lightgreen", high="darkgreen", 
                                     guide="colorbar",na.value="white")+
    theme_classic()
}

# We are not able to extract the value for all sites so we assign the nearest value (TO BE DISCUSSED)  



x <- read_csv("remainingPlots_elev.csv")
plotLocations$test[is.na(plotLocations$test)] <- x$test[match(plotLocations$PlotObservationID[is.na(plotLocations$test)], x$PlotObservationID)]

remainingPlots<- plotLocations[is.na(plotLocations$test),]
# I was not able to do this with foreach so takes about 1 hour to calculate
for(i in 1:nrow(remainingPlots)){
  # Buffer point locations by desired distance
  ptBuff <- st_buffer(remainingPlots[i,], dist = 10000)
  # Crop raster to buffered point
  inbuff <- crop(elev, vect(ptBuff))
  # Convert to points
  dat <- as.points(inbuff)
  dat <- st_as_sf(dat)
  # If all values are NA the dataframe will be empty, move to next
  if(nrow(dat) == 0){
    next
  }
  ptdist <- st_distance(remainingPlots[i,], dat)
  dat$pdistance <- as.numeric(ptdist)
  m <- as.numeric(min(ptdist))
  val <- dat %>% 
    dplyr::filter(pdistance == m) %>% 
    select(Elev) %>% 
    st_drop_geometry()
  val <- as.numeric(val)
  remainingPlots[i,]$test <- val
}

# assign the new values to the previous NAs  
plotLocations[is.na(plotLocations$test),] <- remainingPlots
#write_csv(plotLocations, "remainingPlots_elev.csv")


# Make dataset with all plots that are still not assigned --> probably due to the raster nature --> maybe make the distance larger
remainingPlots<- plotLocations[is.na(plotLocations$test),]
# does it makes sense? Looking at the coordinates it is fair to say that they are at sea level
remainingPlots$test <- "0"
plotLocations[is.na(plotLocations$test),] <- remainingPlots
# Plot these remaining not extracted plots --> all close to the coast
if(plotting){
  ggplot() +   
    geom_spatraster(data = elev) + 
    geom_sf(data = remainingPlots, color = "black", size = 1) + 
    coord_sf()+theme_minimal()+scale_fill_continuous(low="lightgreen", high="darkgreen", 
                                                     guide="colorbar",na.value="white")
}

# prepare for saving
saving= T
if(saving){
  st_geometry(plotLocations) <- NULL
  colnames(plotLocations)[6]<- "elev"
  fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
  fullPlotData<- left_join(fullPlotData, plotLocations)
  fullPlotData$elev <- as.numeric(fullPlotData$elev)
  write_csv(fullPlotData, "fullPlotData_ESy.csv")
}
 

###### 4.4 Remove data ###### 
eva <- read_csv("fullPlotEva_ESy.csv")
eva <- eva[!eva$PlotObservationID %in% fullPlotData$PlotObservationID[is.na(fullPlotData$hfp)],]
fullPlotData <- fullPlotData[!is.na(fullPlotData$hfp),]

# Save
saving= T
if(saving){
  write_csv(fullPlotData, "fullPlotData_ESy.csv")
  write_csv(eva, "fullPlotEva_ESy.csv")
  
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#### 5 HEADER ####
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
###### 5.1 Prepare #####
#------------------------------------------------------------------------------#
# read in ecoregion data
eco<- st_read("../EIVE Data/Ecoregions/Ecoregions2017.shp")
# transform CSR
eco <- st_transform(eco, CRS("+proj=longlat +datum=WGS84"))

# look at data and see which realm Europe belongs to
head(eco)
unique(eco$REALM)

# only get Europe realm
eco <- eco[eco$REALM =="Palearctic",]

# get euro-med regions
medRegions <- read_sf("../Europe-regions-shapefiles-2023", "Emed_regions")
medRegions <- st_transform(medRegions, CRS("+proj=longlat +datum=WGS84"))

# create box around Europe
europe_bbox <- st_bbox(c(xmin = -30, xmax = 80, ymin = 25, ymax = 90))
europe_polygon <- st_as_sfc(europe_bbox)

# crop to regions
sf_use_s2(F)
eco <- st_crop(eco, europe_polygon)
eco <- st_crop(eco, medRegions, crop=F)
sf_use_s2(F)


# plot
ggplot() + 
  geom_sf(data = eco, size = 1.5, aes(fill = COLOR),show.legend = FALSE) + 
  coord_sf()+
  theme_minimal()

# read in fullPlotData
fullPlotData <- read_csv("fullPlotData_ESy.csv")
fullPlotData <- fullPlotData[runif(length(fullPlotData$PlotObservationID))>0,]

# make spatial object
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5)]
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")


sf_use_s2(T) #S2 can provide better performance for certain types of spatial operations, especially when dealing with large spatial datasets.
joinedData <- st_join(plotLocations, eco, join = st_nearest_feature)

joinedData <- joinedData[, c("PlotObservationID","Country","Longitude","Latitude","ECO_NAME","ECO_ID","COLOR", "geometry")]


#------------------------------------------------------------------------------#
###### 5.2 Bohn #####
#------------------------------------------------------------------------------#
bohn <- st_read("../EIVE Data/ESy-master/data/DUNES_BOHN/Dunes_BohnMap_buffer500m.shp")

# transform CSR
bohn <- st_transform(bohn, CRS("+proj=longlat +datum=WGS84"))

# look at data and see which realm Europe belongs to
head(bohn)
unique(bohn$CODE)

# P1-2 Arctic
# P3-4 Baltic
# P5-8 Arctic
# P9-12 Mediteranean
# P13-16 Black Sea

# plot
ggplot() + 
  geom_sf(data = bohn, size = 5, aes(fill = CODE),show.legend = T) + 
  coord_sf()+
  theme_minimal()

sf_use_s2(T) #S2 can provide better performance for certain types of spatial operations, especially when dealing with large spatial datasets.
joinedData <- st_join(joinedData, bohn, join = st_within)

joinedData <- joinedData[, c(1:8, 16,25)]

###### 5.3 Coast #####
coast <- st_read("../EIVE Data/Coast/Europe_coastline.shp")
coast <- st_transform(coast, CRS("+proj=longlat +datum=WGS84"))
head(coast)

# buffer 10 km
coast <- st_buffer(coast, dist= 10000)

# plot
ggplot() + 
  geom_sf(data = coast, size = 1.5,show.legend = FALSE) + 
  coord_sf()+
  theme_minimal()

# create list of all plots that are within this dataset
coast_plots <- st_join(joinedData, coast,join= st_within)
coast_plots$Shape_Leng[!is.na(coast_plots$Shape_Leng)]<- 1


# assign to nearest coast value
coast_test <- st_join(coast_plots[!is.na(coast_plots$Shape_Leng),-c(9:11)], bohn, join= st_nearest_feature)

joinedData$CODE[joinedData$PlotObservationID %in% coast_test$PlotObservationID] <- coast_test$CODE[match(joinedData$PlotObservationID[joinedData$PlotObservationID %in% coast_test$PlotObservationID], coast_test$PlotObservationID)]

joinedData <- joinedData |> mutate(COAST_TYPE = case_when(
  CODE %in% c("P1", "P2") ~ "ARC_COAST",  # P1-P2 -> ARC_COAST
  CODE %in% c("P3", "P4") ~ "BAL_COAST",  # P3-P4 -> BAL_COAST
  CODE %in% c("P5", "P6", "P7", "P8") ~ "ATL_COAST",  # P5-P8 -> ATL_COAST
  CODE %in% c("P9", "P10", "P11", "P12") ~ "MED_COAST",  # P9-P12 -> MED_COAST
  CODE %in% c("P13", "P14", "P15", "P16") ~ "BLA_COAST",  # P13-P16 -> BLA_COAST
  is.na(CODE) ~ "N_COAST",  # NA values -> N_COAST
  TRUE ~ "N_COAST"  # In case there are any unexpected values
))

st_geometry(plotLocations) <- NULL
st_geometry(joinedData) <- NULL
fullPlotData<- left_join(fullPlotData, joinedData)

#write_csv(fullPlotData, "fullPlotData_EUNIS.csv")

#### 6 SAC ####
library(devtools)
library(moranfast)
fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
moranfast(fullPlotData$ENS0, fullPlotData$Longitude, fullPlotData$Latitude)
