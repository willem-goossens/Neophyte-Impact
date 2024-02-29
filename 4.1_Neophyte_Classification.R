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
extra<-setdiff(all, neophyteNamesEU)

###### classify ######
# neophyteDef 
neophyteDefEU$Neophyte[neophyteDefEU$species %in% neophyteNames]<- "extra"
length(unique(neophyteDefEU$species[neophyteDefEU$Neophyte=="extra"]))
# we miss some species (more in neophyte names --> 25 species from not_defined --> are present in old document but not in new) 
# so we test it on the large database
eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% neophyteNames] <- "extra"
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte =="extra"]))
# one extra --> NA


