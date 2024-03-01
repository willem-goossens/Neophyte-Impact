rm(list=ls())

##### 1 LOAD ####
install.packages("tmap", repos = c("https://r-tmap.r-universe.dev", "https://cloud.r-project.org"))
library(tmap)   
library(leaflet) 
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

##### 2 DATA ####
###### 2.1 neophytes ######
neophyte<- readxl::read_excel("../Neophyte-Impact/country_species-2024-01-30-IA-VK.xlsx", sheet="country_speciesWillem")


###### 2.2 fullPlotData #####
# Load and make smaller to conocate
fullPlotEva <- read_csv("fullPlotEva.csv")
eva2<- fullPlotEva[,c("PlotObservationID","species")]
fullPlotData<- read_csv("fullPlotData.csv")
fullPlot2<- fullPlotData[,c("PlotObservationID","Region")]

# Join eva and header (species per plot and plot info respectively)
eva_country<- left_join(eva2, fullPlot2)
# Join eva_country with neophyte data from Irena --> gives all alien species per country
eva_country_neophyte<- left_join(eva_country, neophyte, by= c("Region"= "Region", "species"= "species"))
eva_country_neophyte<- subset(eva_country_neophyte, select=-c(SeqID, FloraVegSpecies))

# Take dataset with all neophytes within EU
neophyteDefEU<- (eva_country_neophyte[(!is.na(eva_country_neophyte$Neophyte)),])
# Retrieve all names of all alien species to Europe
neophyteNamesEU<- unique(neophyteDefEU$species)
# In total, we have 1772 alien species in European countries

# Check whether joining worked fine
length(unique(neophyteDefEU$species[neophyteDefEU$Region=="Poland"]))
length(unique(neophyte$species[neophyte$Region== "Poland"& neophyte$Neophyte=="neo"]))
# Difference --> there is a NA species in the second --> but not if checked --> weird that it comes up --> CHECK!


###### 2.3 extra-EU #######
# Data on which species are neophytes from outside of Europe
neophyteDefinitions <- read_csv("../Neophyte Assignments/UniqueTaxaEurope-2023-04-23.csv")
# Assign these names to the eva list
neophyteNames <- neophyteDefinitions$species[neophyteDefinitions$statusEurope == "neo"]
neophyteNames <- unique(neophyteNames[neophyteNames %in% fullPlotEva$species])

# Some species should be excluded --> Check which 
exclude <- neophyteDefinitions$species[neophyteDefinitions$statusEurope == "exclude"]
exclude <- exclude[exclude %in% fullPlotEva$species]
# Maybe best to remove these from eva --> otherwise bias in application of names
# But first check how they are incorporated in the new file --> not neophyte
neophyteDefEU[(neophyteDefEU$species %in% exclude),]
eva_country_neophyte <- eva_country_neophyte[!(eva_country_neophyte$species %in% exclude),]




#### 3 ANALYSIS ####
###### 3.1 intra-EU ######
# Take names all intra-EU species --> also possible (setdiff(neophyteNamesEU, neophyteNames))
intra_EU<- neophyteNamesEU[!(neophyteNamesEU %in% neophyteNames)]

# The sums do not add up to all species in neophyteNamesEU, we check which species are present in neophyteNames and not in neophyteNamesEU
all<- c(intra_EU, neophyteNames)
not_defined<-setdiff(all, neophyteNamesEU) # we checked other way around, is empty character vector: setdiff(neophyteNamesEU, all)
# CONCLUSION: there are 25 species that were defined extra-European in the old dataset but were now included in the list as native, causing the
# mismatch between both datasets. 

# Extra test to see where they can be found
country_species_number<- eva_country_neophyte |> group_by(Region, species, Neophyte) |> summarise( n= n())
print(country_species_number[country_species_number$species %in% not_defined, ], n= 48)
# most species are found in Turkey --> in old one defined as extra-European but now included.
print(neophyte[neophyte$species %in% not_defined,], n=48)
# Citrus x limon everywhere alien
# Cucumis melo everywhere alien
# Cuscuta epilinum everywhere alien
# Galanthus elwesii alien in Romania (KEW, CABI not)
# Oplismenus hirtellus subsp. undulatifoliu alien in Georgia, Croatia, Switzerland, Slovenia, Turkey (everywhere alien CABI)
# Phoenix dactylifera alien

# TO DO 1: ISLANDS CHECK UP
# TO DO 2: TURKEY AND GEORGIA CHECK

# Remove the species for which we are certain that they are from Turkey (based on KEW (and minor proportion CABI) from NeophyteNames)
# Only for species not present in other countries except for Turkey of course
Turkey_not_neophyte <- c("Camelina laxa","Fritillaria persica","Crocus kotschyanus subsp. suworowianus","Lepyrodiclis holosteoides", 
                         "Moluccella laevis","Potentilla divaricata", "Roemeria refracta","Sisymbrium septulatum",
                         "Sorghum halepense var. muticum", "Sternbergia vernalis","Thlaspi huetii","Trigonella capitata", 
                         "Tripleurospermum decipiens", "Tripleurospermum disciforme", "Vicia noeana"  )
Georgia_not_neophye<- c("Veronica ceratocarpa", "Schoenoplectus juncoides")
# we check our previous number of strange species against these only in countries where they are native
setdiff(not_defined, c(Turkey_not_neophyte,Georgia_not_neophye))
# There is one species that was found to be native to spain according to KEW
others_not_neophyte<- c("Delphinium obcordatum")

# We remove all species that were wrongly indicated in the first dataset from our list
extra_EU<- neophyteNames[!(neophyteNames %in% c(Turkey_not_neophyte, Georgia_not_neophye, others_not_neophyte))]
# These are all neophytes that are extra-European
# We now retrieve the six species that have to be checked
all<- c(intra_EU, neophyteNames)
not_defined<-setdiff(all, neophyteNamesEU)
# For now, we will just exclude these from our dataset
eva_country_neophyte<- eva_country_neophyte[!(eva_country_neophyte$species %in% not_defined),]


###### 3.2 re-classify ######
# We now reclassify the data of eva_country to the newest version, where all species from the old file (extra-EU) are categorised as extra
# The rest of the species are considered intra-EU neophytes
eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% extra_EU] <- "extra"
eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% intra_EU] <- "intra"
eva_country_neophyte$Neophyte[is.na(eva_country_neophyte$Neophyte)] <- "native"
unique(eva_country_neophyte$Neophyte)

# Check number of intra and extra European species and all natives
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte=="extra"]))
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte=="intra"]))
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte=="native"]))

# Summarise again now only number and region in order to see how the data varies in European countries
country_neophyte<- eva_country_neophyte |> distinct(Region, Neophyte, species) %>% group_by(Region, Neophyte) |> summarise(n=n())
table(country_neophyte['Neophyte'])

##### 4 MAP ####
###### 4.1 MED regions #####
# We load the Med regions and give correct CRS
medRegions <- read_sf("../Europe-regions-shapefiles-2023", "Emed_regions")
medRegions <- st_transform(medRegions, CRS("+proj=longlat +datum=WGS84"))


###### 4.2 absolute #####
# We extract the data for only the extra EU species and map this
country_extra<- country_neophyte[country_neophyte$Neophyte=="extra",]
Europe<- left_join(medRegions,country_extra,  by= c("Region"="Region"))
extra_EU<- ggplot()+
  geom_sf(data= Europe, aes(fill= n))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Extra-European Neophyte Distribution in European Regions") +
  theme_minimal()
ggsave(extra_EU, file="extra_EU.png", bg="white")

# We extract the data for only the intra EU species and map this
country_intra<- country_neophyte[country_neophyte$Neophyte=="intra",]
Europe_in<- left_join(medRegions,country_intra,  by= c("Region"="Region"))
intra_EU<- ggplot()+
  geom_sf(data= Europe_in, aes(fill= n))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Intra-European Neophyte Distribution in European Regions") +
  theme_minimal()
ggsave(intra_EU, file="intra_EU.png", bg="white")

# We extract the data for only the native EU species and map this
country_native<- country_neophyte[country_neophyte$Neophyte=="native",]
Europe_nt<- left_join(medRegions,country_native,  by= c("Region"="Region"))
native_EU<-ggplot()+
  geom_sf(data= Europe_nt, aes(fill= n))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Native Distribution in European Regions") +
  theme_minimal()
ggsave(native_EU, file="native_EU.png", bg="white")


###### 4.2 relative #####
# Group by region and calculate the relative proportions of all species
country_neophyte_relative <- eva_country_neophyte %>%
  distinct(Region, Neophyte, species) %>%  
  group_by(Region, Neophyte) %>%
  summarise(n = n_distinct(species)) %>%
  group_by(Region) %>%
  mutate(relative_proportion = n / sum(n))

# We extract the data for only the extra EU species and map this
country_extra_rel<- country_neophyte_relative[country_neophyte_relative$Neophyte=="extra",]
Europe<- left_join(medRegions,country_extra_rel,  by= c("Region"="Region"))
extra_EU_rel<- ggplot()+
  geom_sf(data= Europe, aes(fill= relative_proportion))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Extra-European Neophyte Distribution in European Regions") +
  theme_minimal()
ggsave(extra_EU_rel, file="extra_EU_rel.png", bg="white")

# We extract the data for only the intra EU species and map this
country_intra_rel<- country_neophyte_relative[country_neophyte_relative$Neophyte=="intra",]
Europe_intra<- left_join(medRegions,country_intra_rel,  by= c("Region"="Region"))
intra_EU_rel<- ggplot()+
  geom_sf(data= Europe_intra, aes(fill= relative_proportion))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Intra-European Neophyte Distribution in European Regions") +
  theme_minimal()
ggsave(intra_EU_rel, file="intra_EU_rel.png", bg="white")

# We extract the data for only the native species and map this
country_native_rel<- country_neophyte_relative[country_neophyte_relative$Neophyte=="native",]
Europe_native<- left_join(medRegions,country_native_rel,  by= c("Region"="Region"))
native_EU_rel<- ggplot()+
  geom_sf(data= Europe_native, aes(fill= relative_proportion))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Native Distribution in European Regions") +
  theme_minimal()
ggsave(native_EU_rel, file="native_EU_rel.png", bg="white")
