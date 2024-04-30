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
fullPlotData <- read_csv("fullPlotData_cover_all_layer.csv", show_col_types = FALSE)
fullPlotData<- fullPlotData[,c(1:5, 9:10)]

# reduce size
fast <- F
if(fast) {
  fullPlotData <- fullPlotData[runif(length(fullPlotData$PlotObservationID)) > 0.99,]
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
plotLocations[is.na(plotLocations$test),] <- remainingPlots

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
  colnames(plotLocations)[8]<- "hfp2009"
  fullPlotData <- read_csv("fullPlotData_cover_all_layer.csv", show_col_types = FALSE)
  fullPlotData<- left_join(fullPlotData, plotLocations)
  #write_csv(fullPlotData, "fullPlotData2.csv")
}

###### 3.3 HFP 1993 #####
# Get plot locations and give geometry
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5, 9:10, 54)]
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

# We are not able to extract the value for all sites so we assign the nearest value (TO BE DISCUSSED)  
remainingPlots<- plotLocations[is.na(plotLocations$test),]


# I was not able to do this with foreach so takes about 1 hour to calculate
for(i in 1:nrow(remainingPlots)){
  # Buffer point locations by desired distance
  ptBuff <- st_buffer(remainingPlots[i,], dist = 30000)
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
  # fullPlotData <- read_csv("fullPlotData2.csv", show_col_types = FALSE)
  colnames(plotLocations)[8]<- "hfp1993"
  fullPlotData<- left_join(fullPlotData, plotLocations)
  # write_csv(fullPlotData, "fullPlotData2.csv")
}
# Note that the lengths are not of the same length --> strange?




###### 3.4 Time #####
# Load dataset again if not wanting to do all previous work
# fullPlotData<- read_csv("fullPlotData2.csv", show_col_types = FALSE)

# Check if there are any NAs --> not here but yes later --> are dates 00.00.0000 --> remove later!
anyNA(fullPlotData$Date)
# Make dates
fullPlotData$Date<- as.Date(fullPlotData$Date, format = "%d.%m.%Y")

# Assign extracted data to correct plot
fullPlotData$hfp<- ifelse(fullPlotData$Date > dmy("01012001"), fullPlotData$hfp2009, fullPlotData$hfp1993)

# Save
saving= T
if(saving){
  #write_csv(fullPlotData, "fullPlotData_.csv")
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
chelsaT<- terra::mask(chelsaT, medRegions) 
chelsaP<- terra::mask(chelsaP, medRegions) 
elev<- terra::mask(elev, medRegions) 

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
plotLocations <- plotLocations[, c(1:5, 9:10, 55)]
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


# I was not able to do this with foreach so takes about 1 hour to calculate
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
  # fullPlotData <- read_csv("fullPlotData2.csv", show_col_types = FALSE)
  colnames(plotLocations)[8]<- "chelsaT"
  fullPlotData<- left_join(fullPlotData, plotLocations)
  write_csv(fullPlotData, "fullPlotData_cover_all_layer.csv")
}



###### 4.3 P mean #####
# Get plot locations and give geometry
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5, 9:10, 57)]
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
  colnames(plotLocations)[8]<- "chelsaP"
  fullPlotData<- left_join(fullPlotData, plotLocations)
  # write_csv(fullPlotData, "fullPlotData2.csv")
}

###### 4.3 Elev mean #####
# Get plot locations and give geometry
plotLocations <- st_as_sf(fullPlotData, coords = c("Longitude","Latitude"), remove = FALSE)
plotLocations <- plotLocations[, c(1:5, 9:10, 58)]
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

# Make dataset with all plots that are still not assigned --> probably due to the raster nature --> maybe make the distance larger
remainingPlots<- plotLocations[is.na(plotLocations$test),]
# Plot these remaining not extracted plots --> all close to the coast
if(plotting){
  ggplot() +   
    geom_spatraster(data = elev) + 
    geom_sf(data = remainingPlots, color = "black", size = 1) + 
    coord_sf()+theme_minimal()+scale_fill_continuous(low="lightgreen", high="darkgreen", 
                                                     guide="colorbar",na.value="white")
}

# prepare for saving
saving= F
if(saving){
  st_geometry(plotLocations) <- NULL
  # fullPlotData <- read_csv("fullPlotData2.csv", show_col_types = FALSE)
  colnames(plotLocations)[8]<- "elev"
  fullPlotData<- left_join(fullPlotData, plotLocations)
  #write_csv(fullPlotData, "fullPlotData_cover_all_layer.csv")
}


fullPlotDataNas <- is.na(fullPlotData)
# Get the plot observations
plotsForWhichNotAllIndicatorValuesWereEstimated <- fullPlotDataNas[,1]
# 53
for(i in 53:dim(fullPlotDataNas)[2]) {
  plotsForWhichNotAllIndicatorValuesWereEstimated <- plotsForWhichNotAllIndicatorValuesWereEstimated | fullPlotDataNas[,i]
}
numberOfPlotsWithMissingIndicatorValues <- sum(plotsForWhichNotAllIndicatorValuesWereEstimated)

fullPlotData<- fullPlotData[!plotsForWhichNotAllIndicatorValuesWereEstimated,]
#write_csv(fullPlotData, "fullPlotData_cover_all_layer_cleaned.csv")

eva <- read_csv("fullPlotEva_cover_all_layer.csv")
eva<- eva[eva$PlotObservationID %in% fullPlotData$PlotObservationID,]
#write_csv(eva,"fullPlotEva_cover_all_layer_cleaned.csv")
