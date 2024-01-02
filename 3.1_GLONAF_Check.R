rm(list=ls())

# Load lirbrary
library(dplyr)
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

# Load data Zang et al 2023
native<-read.csv("../Neophyte Assignments/wf_dat.csv")
neophyte<- read.csv("../Neophyte Assignments/glonaf_dat.csv")

# Load our cleaned up header data
fullPlotData<- read.csv("fullPlotData.csv")

# Check how data is organised and retrieve only data on Europe
# We also checked for Australasia, but this is only Australia, hence, we can only keep the Europe part of the dataset for our analysis
head(native)
unique(native$botanical_continent)
EUnative<- native[native$botanical_continent=="EUROPE",]

# Same for neophytes
head(neophyte)
unique(neophyte$botanical_continent)
EUneo<- neophyte[neophyte$botanical_continent=="EUROPE",]

# Check the countries of the Zang data against our data
WF_regions<-unique(EUnative$geo_entity)
MED_regions<-unique(fullPlotData$Region)
GLONAF_regions<- unique(EUneo$geo_entity)

# Correct countries for the two datasets
correctCountries<- data.frame(Med=c("Rf.NW", "Rf.N","Rf.E", "Rf.C","Rf.S", "Luxemburg","Bosnia.Herzegovina", "Italy" ,"Czech.Republic", "Greece", "France", "Spain", 
                                    "Portugal" , "Moldavia"), 
                              WF=c("Northwest European Russia", "North European Russia", "East European Russia", "Central European Russia", "South European Russia", 
                                "Luxembourg", "Bosnia and Herzegovina", "Italy excl. Sardinia and Sicily", "Czech Republic", "Greece excl. Crete and East Aegean",
                                "France incl. Channel Islands and Monaco excl. Corse","Spain mainland", "Portugal mainland", "Moldova"))

# Correct WF dataset
index<- WF_regions %in% correctCountries$WF
WF_regions[index]<- correctCountries$Med[match(WF_regions[index], correctCountries$WF)]

# Correct GLONAD dataset
index<- GLONAF_regions %in% correctCountries$WF
GLONAF_regions[index]<- correctCountries$Med[match(GLONAF_regions[index], correctCountries$WF)]

# Check which are still different WF
intersect(WF_regions, MED_regions)
print(sort(setdiff(WF_regions, MED_regions)))
print(setdiff( MED_regions, WF_regions))

# Check which are the same
intersect(GLONAF_regions, MED_regions)
print(sort(setdiff(GLONAF_regions, MED_regions)))
print(setdiff( MED_regions, GLONAF_regions))


# Zang et al excluded the islands
# Load all data
allPlotsWithRegion<- read.csv("plot_to_region_assignment.csv")
# Load header data
header <- read_delim("../EVA Data/171_NeophyteInvasions20230216_notJUICE_header.csv", "\t")
header <- header[!(is.na(header$Latitude) | is.na(header$Longitude)),]
plotLocations <- st_as_sf(header, coords = c("Longitude","Latitude"), remove = FALSE)
st_crs(plotLocations) <- CRS("+proj=longlat +datum=WGS84")
# Load medRegions
medRegions <- read_sf("../Europe-regions-shapefiles-2023", "Emed_regions")
medRegions <- st_transform(medRegions, CRS("+proj=longlat +datum=WGS84"))

# Make the boxes for each country, as well as for Europe
boundingBoxes <- medRegions
for(i in 1:length(boundingBoxes$country)) {
  boundingBoxes$geometry[i] <- st_as_sfc(st_bbox(boundingBoxes$geometry[i]))
}
EUbox<-st_as_sfc(st_bbox(boundingBoxes$geometry))

# Code for plotting the plots and regions
regionToLookAt <- ("Serbia+Kosovo") # With this line you can select which region to plot
regionShape<- medRegions[medRegions$Region==regionToLookAt,]
boundingBox <- boundingBoxes$geometry[match(regionToLookAt, medRegions$Region)]
plotsInRegion <- plotLocations[allPlotsWithRegion$Region == regionToLookAt & runif(length(header$Country)) > 0,] # Here you can reduce the number of plots shown

# Make sf
medRegions.df <- fortify(medRegions)
regionShape <- fortify(regionShape)
boundingBox <- fortify(boundingBox)

# Plot
ggplot() +   
  #geom_sf(data = boundingBox, color = "red") + 
  geom_sf(data = regionShape, color = "black", fill = "cyan1") + 
  geom_sf(data = plotsInRegion, color = "red", size = 0.5) + 
  coord_sf()


# Glonaf
glonafRegions <- read_sf("../GloNAF_Shapefile", "regions2")
glonafRegions <- st_transform(glonafRegions, CRS("+proj=longlat +datum=WGS84"))

# Whole world
glonafRegions<- fortify(glonafRegions)
ggplot()+geom_sf(data=EUbox, color="red")+
  geom_sf(data=glonafRegions, color="black", fill="cyan")

# Only Europe
glonafRegionList<- read.csv("../GloNAF_Shapefile/Region_GloNAF_vanKleunenetal2018Ecology.csv")
glonafRegionList<- glonafRegionList[glonafRegionList$tdwg1_name=="Europe", ]  
unique(glonafRegionList$tdwg4_name)
index<- glonafRegions$OBJIDsic %in% glonafRegionList$OBJIDsic
newGlonaf<- glonafRegions[index,]
newGlonaf<- left_join(newGlonaf, glonafRegionList, by= "OBJIDsic")
newGlonaf<- fortify(newGlonaf)

# Plot European Glonaf plots
# Europe_Plots<- fortify(plotLocations)
ggplot()+geom_sf(data=EUbox, color="red")+
  geom_sf(data=newGlonaf, color="black", fill="cyan")+
 # geom_sf(data = Europe_Plots, color = "red", size = 0.5) + 
   coord_sf()