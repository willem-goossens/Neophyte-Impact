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
###### neophytes ######
neophyte<- readxl::read_excel("../Neophyte-Impact/country_species-2024-01-30-IA-VK.xlsx", sheet="country_speciesWillem")

###### fullPlotData #####
# Load and make smaller to conocate
fullPlotEva <- read_csv("fullPlotEva.csv")
eva2<- fullPlotEva[,c("PlotObservationID","species")]
fullPlotData<- read_csv("fullPlotData.csv")
fullPlot2<- fullPlotData[,c("PlotObservationID","Region")]

# Join
eva_country<- left_join(eva2, fullPlot2)
eva_country_neophyte<- left_join(eva_country, neophyte, by= c("Region"= "Region", "species"= "species"))

# Take dataset with all neophytes within EU
neophyteDefEU<- (eva_country_neophyte[(!is.na(eva_country_neophyte$Neophyte)),])
neophyteNamesEU<- unique(neophyteDefEU$species)

# Check whether joining worked fine
length(unique(neophyteDefEU$species[neophyteDefEU$Region=="Poland"]))
length(unique(neophyte$species[neophyte$Region== "Poland"& neophyte$Neophyte=="neo"]))

# Data on which species are neophytes from outside of Europe
neophyteDefinitions <- read_csv("../Neophyte Assignments/UniqueTaxaEurope-2023-04-23.csv")
# Assign these names to the eva list
neophyteNames <- neophyteDefinitions$species[neophyteDefinitions$statusEurope == "neo"]
exclude <- neophyteDefinitions$species[neophyteDefinitions$statusEurope == "exclude"]
exclude <- exclude[exclude %in% fullPlotEva$species]

neophyteNames <- neophyteNames[neophyteNames %in% fullPlotEva$species]
# Take names all intra-EU species 
# Also possible (setdiff(neophyteNamesEU, neophyteNames))
intra_EU<- neophyteNamesEU[!(neophyteNamesEU %in% neophyteNames)]

# The sums do not add up to all species in neophyteNamesEU, we check which species are present in neophyteNames and not in neophyteNamesEU
all<- c(intra_EU, neophyteNames)
not_defined<-setdiff(all, neophyteNamesEU) # we checked other way around, is empty character vector: setdiff(neophyteNamesEU, all)

# There are duplicates in neophyteNames, we remove them and check again whether the sum checks out now
neophyteNames<- unique(neophyteNames)
all<- c(intra_EU, neophyteNames)
not_defined <-setdiff(all, neophyteNamesEU)

###### classify ######
# neophyteDef --> only species with no NA --> but too little --> some species are present in old file of Irena and not in new
neophyteDefEU$Neophyte[neophyteDefEU$species %in% neophyteNames]<- "extra"
length(unique(neophyteDefEU$species[neophyteDefEU$Neophyte=="extra"]))
# we miss some species (more in neophyte names --> 25 species from not_defined --> are present in old document but not in new) 
# so we test it on the large database
eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% neophyteNames] <- "extra"
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte =="extra"]))
# one extra --> NA

# We look at old species in not_defined file --> where are they present
writeClipboard(not_defined)
# We checked species manually
#   Camelina laxa --> Turkey only
#   Citrus x limon --> no native range, introduced in Spain, Portugal and Albania
#   Cucumis melo --> alien and introduced
#   Cuscuta epilinum --> Turkey
#   Fritillaria persica --> Turkey
#   Lepyrodiclis holosteoides --> Turkey
#   Tripleurospermum disciforme --> Turkey
#   Moluccella laevis --> Turkey
#   Oplismenus hirtellus subsp. undulatifolius --> alien --> but present 
#   Phoenix dactylifera --> introduced to Turkey but not native
#   Potentilla divaricata --> Turkey
#   Roemeria refracta --> Turkey
#   Schoenoplectus juncoides --> alien --> introduced Portugal
#   Sisymbrium septulatum -_> Turkey
#   Sorghum halepense var. muticum --> Turkey
#   Sternbergia vernalis --> Turkey
#   Trigonella capitata --> Turkey
#   Tripleurospermum decipiens --> Turkey
#   Veronica ceratocarpa --> Turkey
#   Vicia noeana --> Turkey
#   Crocus kotschyanus subsp. suworowianus --> Turkey
#   Galanthus elwesii --> Albania, Bulgaria, Greece, Ukraine, Yugoslavia, Turkey
#   Thlaspi huetii --> Turkey

# Majority is only from turkey --> check in eva_country_neophyt
Turkey<- eva_country_neophyte[eva_country_neophyte$Region=="Turkey",]
unique(Turkey$species[Turkey$Neophyte=="neo"]) 

# Check which species are also found in Turkey
Turkey_not_defined<- c()
for (i in 1:25){
y <- (which(eva_country_neophyte$species[eva_country_neophyte$Region=="Turkey"]== not_defined[i]))
  if (length(y)>0) {
    Turkey_not_defined[i]<- not_defined[i]
  }
}

print(Turkey_not_defined)
# Look at how they are defined now --> extra -->incorrect --> should be native for most of them --> check for which
Turkey$Neophyte[(Turkey$species %in% Turkey_not_defined)]
# But we need to be sure that they are not found in other countries
rest_not_defined<- c()
for (i in 1:25){
  y <- (which(eva_country_neophyte$species[!eva_country_neophyte$Region=="Turkey"]== not_defined[i]))
  if (length(y)>0) {
    rest_not_defined[i]<- not_defined[i]
  }
}  
print(rest_not_defined)
both_not_defined<- setdiff(rest_not_defined, Turkey_not_defined)
print(both_not_defined)


# Extra test to see where they can be found
test<- eva_country_neophyte |> group_by(Region, species, Neophyte) |> summarise( n= n())
print(test[test$species %in% not_defined, ], n= 48)
# 
print(neophyte[neophyte$species %in% not_defined,], n=48)
# Citrus x limon everywhere alien
# Cucumis melo everywhere alien
# Cuscuta epilinum everywhere alien
# Galanthus elwesii alien in Romania (KEW, CABI not)
# Oplismenus hirtellus subsp. undulatifoliu alien in Georgia, Croatia, Switzerland, Slovenia, Turkey (everywhere alien CABI)
# Phoenix dactylifera alien

country<- test|> group_by(Region, Neophyte) |> summarise(n=n())

medRegions <- read_sf("../Europe-regions-shapefiles-2023", "Emed_regions")
medRegions <- st_transform(medRegions, CRS("+proj=longlat +datum=WGS84"))

Europe<- left_join(country, medRegions, by= c("Region"="Region"))
Europe<- st_combine(Europe)

Europe_GG<- fortify(Europe)
ggplot()+
  geom_sf(data=Europe_GG)

allPlotsWithRegion[(is.na(allPlotsWithRegion$Region)),]


world_map <- map_data("world")

Europe<- st_combine(medRegions)
Europe_Plots<- fortify(plotLocations)

Europe_GG<- fortify(Europe)
ggplot()+
  geom_sf(data=Europe_GG)+
  geom_sf(data= Europe_Plots, color="red", size=0.5)

allPlotsWithRegion[(is.na(allPlotsWithRegion$Region)),]# List of European countries
european_countries <- c(
  "Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herzegovina",
  "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia",
  "Finland", "France", "Germany", "Gibraltar", "Greece", "Hungary", "Iceland",
  "Ireland", "Italy", "Kosovo", "Latvia", "Liechtenstein", "Lithuania", "Luxembourg",
  "Macedonia", "Malta", "Moldova", "Monaco", "Montenegro", "Netherlands",  "Poland","Norway",
  "Portugal", "Romania", "San Marino", "Serbia", "Slovakia", "Slovenia", "Sweden","Spain", "Switzerland", "Ukraine", "United Kingdom", "Turkey"
)

# Subset the data for European countries
europe_map <- subset(world_map, region %in% european_countries)

# Highlight France and Spain
highlighted_countries <- c("France", "Spain")

# Create a basic map
ggplot(europe_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightgrey", color = "black") +
  geom_polygon(data = subset(europe_map, region %in% highlighted_countries),
               fill = "darkgrey", color = "black") +
  theme_void() 
  
