rm(list=ls())

##### 1 LOAD ####
library(tmap)   
library(leaflet) 
library(readr)
library(ggplot2)
library(sf)
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
neophyte <- read_csv("neophyte_euro.csv", show_col_types = FALSE)


###### 2.2 Eva and Header #####
# Load and make smaller to conocate
fullPlotEva <- read_csv("fullPlotEva_euro.csv", show_col_types = FALSE)
eva2<- fullPlotEva[,c("PlotObservationID","species", "irena","Matched concept")]
fullPlotData<- read_csv("fullPlotData_euro.csv", show_col_types = FALSE)
fullPlot2<- fullPlotData[,c("PlotObservationID","Region")]

# Load summarized data
species_country <- read_csv("country_species_euro_cover.csv", show_col_types = FALSE)

# Join eva and header (species per plot and plot info respectively)
eva_country<- left_join(eva2, fullPlot2)
# Join eva_country with neophyte data from Irena --> gives all alien species per country
eva_country_neophyte<- left_join(eva_country, species_country, by= c("Region"= "Region", "species"= "species", "irena"="irena","Matched concept"= "Matched concept"))
eva_country_neophyte<- subset(eva_country_neophyte, select=-c(n))
head(eva_country_neophyte)

# Take dataset with all neophytes within EU
neophyteDefEU<- (eva_country_neophyte[(!(eva_country_neophyte$Neophyte=="native")& !is.na(eva_country_neophyte$Neophyte)),])
# Retrieve all names of all alien species to Europe
neophyteNamesEU<- unique(neophyteDefEU$species)
#In total, we have 1772 alien species in European countries (before removing NA dates, after doing so 1686, (1699))

# Check whether joining worked fine
Poland <- unique(sort(neophyteDefEU$species[neophyteDefEU$Region=="Poland"]))
Poland2 <-(unique(neophyte$species[neophyte$Region== "Poland"& neophyte$Neophyte=="neo"]))
setdiff(Poland, Poland2)
setdiff(Poland2, Poland)

# Difference --> there is a NA species in the second --> but not if checked --> weird that it comes up --> CHECK!


###### 2.3 extra-EU #######
# Data on which species are neophytes from outside of Europe
neophyteDefinitions <- read_csv("../Neophyte-Impact/Neophyte Assignments/UniqueTaxaEurope-2023-04-23.csv", show_col_types = FALSE)

# get names eva
eva_names <- unique(eva2[, c(2:4)])

# change names to all accepted
neophyteDefinitions$name <- eva_names$species[match(neophyteDefinitions$Matched.concept, eva_names$species)]
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$Matched.concept[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$irena)]
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$Matched.concept[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$`Matched concept`)]
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$species[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$species)]
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$species[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$irena)]
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$species[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$`Matched concept`)]


neophyteDefinitions <- neophyteDefinitions[, c(7, 4,6,3)]
colnames(neophyteDefinitions)<- c("species","name","neophyte","exclude")
neophyteDefinitions <- neophyteDefinitions[!duplicated(neophyteDefinitions),]
# check duplicates

# we check again whether our change was successful
x <- neophyteDefinitions[(duplicated(neophyteDefinitions[,c(1)])| duplicated(neophyteDefinitions[,c(1)], fromLast=TRUE)),]
x <- x[!(duplicated(x[,c(1,3)]) | duplicated(x[,c(1,3)], fromLast=TRUE)),]

change<- data.frame(species= c("Avena sativa","Cynara scolymus", "Fragaria moschata","Gnaphalium",
                               "Raphanus raphanistrum subsp. raphanistrum","Rubus canadensis"),
                    status= c("arch","native","native","native","native","neo"))

# assign to neophyte
neophyteDefinitions$neophyte[neophyteDefinitions$species %in% change$species] <- 
  change$status[match(neophyteDefinitions$species[neophyteDefinitions$species %in% change$species],change$species)]


# we check again whether our change was successful
x <- neophyteDefinitions[(duplicated(neophyteDefinitions[,c(1)])| duplicated(neophyteDefinitions[,c(1)], fromLast=TRUE)),]
x <- x[!(duplicated(x[,c(1,3)]) | duplicated(x[,c(1,3)], fromLast=TRUE)),]

# remove all complete duplicates
neophyteDefinitions <- neophyteDefinitions[!(duplicated(neophyteDefinitions)),]

# check duplicates of species definitions
dup<- neophyteDefinitions[duplicated(neophyteDefinitions$species)| duplicated(neophyteDefinitions, fromLast=TRUE),]
# remove
neophyteDefinitions <- neophyteDefinitions[!duplicated(neophyteDefinitions$species),]
neophyteDefinitions<- neophyteDefinitions[!is.na(neophyteDefinitions$species),]


# join eva names and neophyte 
species <- left_join(eva_names, neophyteDefinitions, by= c("species"="species"))
species[is.na(species$neophyte),] <- left_join(species[is.na(species$neophyte),-c(4:6)], neophyteDefinitions, by= c("irena"="species"))

# 196 --> same as previously
unique(species$species[is.na(species$neophyte)])

# Some species should be excluded --> Check which 
exclude <- neophyteDefinitions$species[neophyteDefinitions$neophyte == "exclude"]
exclude <- exclude[exclude %in% fullPlotEva$species]
# Maybe best to remove these from eva --> otherwise bias in application of names
# But first check how they are incorporated in the new file --> not neophyte
neophyteDefEU <- neophyteDefEU[!(neophyteDefEU$species %in% exclude),]
neophyteNamesEU <- neophyteNamesEU[!(neophyteNamesEU %in% exclude)]
eva_country_neophyte <- eva_country_neophyte[!(eva_country_neophyte$species %in% exclude),]
neophyteDefinitions <- neophyteDefinitions[!(neophyteDefinitions$species %in% exclude),]
neophyteNames <- neophyteDefinitions$species[neophyteDefinitions$neophyte=="neo"]

# Some species are archeophytes
arch <- neophyteDefinitions$species[neophyteDefinitions$neophyte == "arch"]
arch <- arch[arch %in% fullPlotEva$species]
# We will not do anything with them now, but it is possible to adjust the code in the end to also check archeophytes and whether they different (85, 70)


###### 2.4 SUMMARY #####
# All extra EU neophytes defined neo
neophyteDefinitions
neophyteNames
# All aliens and the region they are present
neophyte
neophyteNamesEU


#### 3 ANALYSIS ####
###### 3.1 intra-EU ######
# Take names all intra-EU species --> also possible (setdiff(neophyteNamesEU, neophyteNames)) (899)
intra_EU<- neophyteNamesEU[!(neophyteNamesEU %in% neophyteNames)]

# The sums do not add up to all species in neophyteNamesEU, we check which species are present in neophyteNames and not in neophyteNamesEU (1723)
all<- c(intra_EU, neophyteNames)
not_defined<-setdiff(all, neophyteNamesEU) # we checked other way around, is empty character vector: setdiff(neophyteNamesEU, all)
# CONCLUSION: there are 25 (21, 24) species that were defined extra-European in the old dataset but were now included in the list as native, causing the
# mismatch between both datasets. 

# Extra test to see where they can be found
country_species_number<- eva_country_neophyte |> group_by(Region, species, Neophyte) |> summarise( n= n())
print(country_species_number[country_species_number$species %in% not_defined, ], n= 44)
# most species are found in Turkey --> in old one defined as extra-European but now included.
test<- (neophyte[neophyte$species %in% not_defined,])
test<- test[!(test$species=="Cephalophysis species"),]
# write.csv(test, 'species_not_defined.csv', row.names=FALSE)
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
Turkey_not_neophyte <- c("Camelina laxa","Fritillaria persica","Crocus suworowianus","Lepyrodiclis holosteoides", 
                         "Moluccella laevis","Potentilla divaricata", "Roemeria refracta","Sisymbrium septulatum",
                         "Sorghum halepense var. muticum", "Sternbergia vernalis","Thlaspi huetii","Trigonella capitata", 
                         "Tripleurospermum decipiens", "Tripleurospermum disciforme", "Vicia noeana" , "Gypsophila pilosa","Vicia johannis" )
Georgia_not_neophye<- c("Veronica ceratocarpa", "Schoenoplectus juncoides")
# we check our previous number of strange species against these only in countries where they are native
setdiff(not_defined, c(Turkey_not_neophyte,Georgia_not_neophye))
# There is one species that was found to be native to spain according to KEW
others_not_neophyte<- c("Delphinium obcordatum")
# After checking extra species that can be removed
others_not_neophyte<- c(others_not_neophyte, "Galanthus elwesii", "Galanthus elwesii subsp. elwesii",
                        "Oplismenus hirtellus subsp. undulatifolius")

add_extra<- as.character(c("Citrus x limon","Citrus × limon","Cucumis melo","Phoenix dactylifera"))
add_intra<- as.character(c("Cuscuta epilinum"))

# We remove all species that were wrongly indicated in the first dataset from our list
extra_EU<- neophyteNames[!(neophyteNames %in% c(Turkey_not_neophyte, Georgia_not_neophye, others_not_neophyte))]
extra_EU<- c(extra_EU,add_extra )
# These are all neophytes that are extra-European
# We also add intra-European neophytes
intra_EU<- c(intra_EU, add_intra)

# We now retrieve the six species that have to be checked
all<- c(intra_EU, extra_EU)
not_defined<-setdiff(all,neophyteNamesEU)
not_defined<- c(not_defined, "Cephalophysis species")
not_defined<- setdiff(not_defined, c(add_intra, add_extra))
not_defined


# For now, we will just exclude these from our dataset
eva_country_neophyte<- eva_country_neophyte[!(eva_country_neophyte$species %in% not_defined),]


###### 3.2 re-classify ######
# We now reclassify the data of eva_country to the newest version, where all species from the old file (extra-EU) are categorised as extra
# The rest of the species are considered intra-EU neophytes
#eva_country_neophyte$Neophyte[is.na(eva_country_neophyte$Neophyte)] <- "native"
eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% extra_EU] <- "extra"
eva_country_neophyte$Neophyte[eva_country_neophyte$Neophyte=="neo"]<- "intra"

# unknowns
length(unique(eva_country_neophyte$species[is.na(eva_country_neophyte$Neophyte)]))
unknown <- (eva_country_neophyte[is.na(eva_country_neophyte$Neophyte),])
unknown <- unknown[!duplicated(unknown[, c(2:8)]),]
unknown$species[!unknown$species %in% neophyte$species]

# Check archeophytes
unique(eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% arch])
unique(eva_country_neophyte$species[(eva_country_neophyte$species %in% arch & eva_country_neophyte$Neophyte=="intra")])

# make it native
eva_country_neophyte$Neophyte[eva_country_neophyte$species=="x_Triticosecale rimpaui"]<- "native"

# Check number of intra and extra European species and all natives
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte=="extra"])) #807 (815)
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte=="intra"])) #883 (915 previously) (889)
length(unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte=="native"])) #19327 #19292 (after removing some additional plots) (19617)

# Summarise again now only number and region in order to see how the data varies in European countries
country_neophyte<- eva_country_neophyte |> distinct(Region, Neophyte, species) %>% group_by(Region, Neophyte) |> summarise(n=n())
table(country_neophyte['Neophyte'])

# summarise again based on region and species -_> this time we have extra and intra instead of neo
country_status <- eva_country_neophyte |> group_by(Region, species, irena, `Matched concept`, Neophyte) |> summarise(n=n())

country_status[is.na(country_status$Neophyte),]

# check against old 
#old<- read_csv("species_country_status_new.csv", show_col_types = FALSE)

# check difference with previously made data --> we removed one species (see exclude)
setdiff(country_species_number[, c(1:2,4)], country_status[, c(1:2,4)] )

# Check whether some intra_EU species are native in other regions
native_names<-unique(eva_country_neophyte$species[eva_country_neophyte$Neophyte=="native"])
sum(native_names %in% intra_EU) #818 (before around 848) (823)
native_intra<- native_names[native_names %in% intra_EU]

intra_analysis=F
if(intra_analysis){
# Here we make a new dataframe with those species that are intra in a region in europe as native_intra
# to check whether the effect is just species dependent
eva2<- eva_country_neophyte
eva2$Neophyte[eva2$Neophyte=="native" & eva2$species %in% native_intra]<-"native_intra"
length(unique(eva2$species[eva2$Neophyte=="extra"]))
length(unique(eva2$species[eva2$Neophyte=="intra"]))
length(unique(eva2$species[eva2$Neophyte=="native"])) # 19327-818 = 18509 #18474+818 = 19292
length(unique(eva2$species[eva2$Neophyte=="native_intra"]))

# check whether the classification was successful for all species
sum(unique(eva2$species[eva2$Neophyte=="native"]) %in% intra_EU)


country_neophyte <- eva2 |> distinct(Region, Neophyte, species) %>% group_by(Region, Neophyte) |> summarise(n=n())
table(country_neophyte['Neophyte'])
eva_country_neophyte<- eva2
}

###### 3.3 native SR #####
aggregatedEVA <- eva_country_neophyte |>  group_by(PlotObservationID, Neophyte) |>  summarise(numberOfVascularPlantSpecies = n())

# Only native SR is relevant here (does not make sense to look at the influence of alien species on alien species)
nativeSpR <- aggregatedEVA[aggregatedEVA$Neophyte=="native" | aggregatedEVA$Neophyte=="native_intra" ,]
nativeSpR <- subset(nativeSpR, select= -c(Neophyte))
names(nativeSpR)[names(nativeSpR)=="numberOfVascularPlantSpecies"]<- "nativeSR"

# Are there some plots that do not have any native SR (544) (548)
fullPlotData[!(fullPlotData$PlotObservationID %in% nativeSpR$PlotObservationID),]

# Count number of plots with alien species
plots_with_intra <- sum(aggregatedEVA$Neophyte=="intra", na.rm=T)
plots_with_extra <- sum(aggregatedEVA$Neophyte=="extra", na.rm=T)
plots_with_native<-  sum(aggregatedEVA$Neophyte=="native", na.rm=T)

nativeSR_calculation=F
if(nativeSR_calculation){
# Join native SR with data
fullPlotData<- left_join(fullPlotData, nativeSpR, by="PlotObservationID")
fullPlotData <- fullPlotData |> relocate(nativeSR, .after = numberOfVascularPlantSpecies)
fullPlotData$nativeSR[is.na(fullPlotData$nativeSR)] <- 0
# Check if overview is correct --> difference 19
sum(fullPlotData$numberOfVascularPlantSpecies-fullPlotData$nativeSR)
}

sum(eva_country_neophyte$Neophyte=="extra", na.rm=T)+sum(eva_country_neophyte$Neophyte=="intra", na.rm=T)
# check with natives as well
sum(fullPlotData$numberOfVascularPlantSpecies, na.rm=T)
sum(eva_country_neophyte$Neophyte=="native", na.rm=T)+ sum(eva_country_neophyte$Neophyte=="extra"| eva_country_neophyte$Neophyte=="intra", na.rm=T)

# there is a variance in number of intra and extra species vs the total number of species --> but no NA values
any(is.na(eva_country_neophyte$Neophyte)) 
# difference is also present in difference eva_country and eva_country_neophyte --> check!
# difference is caused by the removal of species with the exclude and not defined part --> check where it is present and reduce with 1
remove<- c(not_defined, exclude)
if(nativeSR_calculation){
  remove_observations <- eva_country$PlotObservationID[eva_country$species %in% remove]
  fullPlotData$numberOfVascularPlantSpecies[fullPlotData$PlotObservationID %in% remove_observations]<-
  fullPlotData$numberOfVascularPlantSpecies[fullPlotData$PlotObservationID %in% remove_observations]-1
}
# checked --> correct

###### 3.4 Save #####
if(intra_analysis){
eva2_country_status<- eva2 |> group_by(Region, species, Neophyte) |> summarise(n=n())
eva2_country_status<- eva2_country_status[,-4]
#write.csv(eva2_country_status,"eva2_country_status_new.csv", row.names = FALSE)
} else {
  species_country_status<- eva_country_neophyte |> group_by(Region, species, Neophyte) |> summarise(n=n())
  species_country_status<- species_country_status[,-4]
  #write.csv(species_country_status,"species_country_status_new.csv", row.names = FALSE)
}

#write.csv(fullPlotData, "fullPlotData_cover_all_layer.csv", row.names=F)

remove<- c(not_defined, exclude)
#write.csv(remove, "not_defined.csv", row.names=FALSE)

#write_csv(new_names, "new_names.csv")
#write_csv(country_status, "species_status_region.csv")

##### 4 MAP ####
###### 4.1 MED regions #####
# We load the Med regions and give correct CRS
medRegions <- read_sf("../Europe-regions-shapefiles-2023", "Emed_regions")
medRegions <- st_transform(medRegions, CRS("+proj=longlat +datum=WGS84"))

density <- left_join(medRegions, country_neophyte, by= c("Region"="Region"))
density$dens<- density$n/density$Shape_Area

density$dens

###### 4.2 absolute #####
# We extract the data for only the extra EU species and map this
country_extra<- country_neophyte[country_neophyte$Neophyte=="extra",]
Europe<- left_join(medRegions,country_extra,  by= c("Region"="Region"))
extra_EU_plot<- ggplot()+
  geom_sf(data= Europe, aes(fill= n))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Extra-European Neophyte Distribution in European Regions") +
  theme_minimal()
#ggsave(extra_EU_plot, file="extra_EU.png", bg="white")

# We extract the data for only the intra EU species and map this
country_intra<- country_neophyte2[country_neophyte$Neophyte=="intra",]
Europe_in<- left_join(medRegions,country_intra,  by= c("Region"="Region"))
intra_EU_plot<- ggplot()+
  geom_sf(data= Europe_in, aes(fill= n))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Intra-European Neophyte Distribution in European Regions") +
  theme_minimal()
#ggsave(intra_EU_plot, file="intra_EU.png", bg="white")

# We extract the data for only the native EU species and map this
country_native<- country_neophyte[country_neophyte$Neophyte=="native",]
Europe_nt<- left_join(medRegions,country_native,  by= c("Region"="Region"))
native_EU<-ggplot()+
  geom_sf(data= Europe_nt, aes(fill= n))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Native Distribution in European Regions") +
  theme_minimal()
#ggsave(native_EU, file="native_EU.png", bg="white")

if(intra_analysis){
  country_native_intra<- country_neophyte[country_neophyte$Neophyte=="native_intra",]
  Europe_nt_in<- left_join(medRegions,country_native_intra,  by= c("Region"="Region"))
  native_intra_EU<-ggplot()+
    geom_sf(data= Europe_nt_in, aes(fill= n))+
    scale_fill_viridis_c(option = "magma",begin = 0.1)+
    labs(title = "Native Distribution in European Regions alien elsewhere") +
    theme_minimal()
  #ggsave(native_intra_EU, file="native_intra_EU.png", bg="white")
}

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
#ggsave(extra_EU_rel, file="extra_EU_rel.png", bg="white")

# We extract the data for only the intra EU species and map this
country_intra_rel<- country_neophyte_relative[country_neophyte_relative$Neophyte=="intra",]
Europe_intra<- left_join(medRegions,country_intra_rel,  by= c("Region"="Region"))
intra_EU_rel<- ggplot()+
  geom_sf(data= Europe_intra, aes(fill= relative_proportion))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Intra-European Neophyte Distribution in European Regions") +
  theme_minimal()
#ggsave(intra_EU_rel, file="intra_EU_rel.png", bg="white")

# We extract the data for only the native species and map this
country_native_rel<- country_neophyte_relative[country_neophyte_relative$Neophyte=="native",]
Europe_native<- left_join(medRegions,country_native_rel,  by= c("Region"="Region"))
native_EU_rel<- ggplot()+
  geom_sf(data= Europe_native, aes(fill= relative_proportion))+
  scale_fill_viridis_c(option = "magma",begin = 0.1)+
  labs(title = "Native Distribution in European Regions") +
  theme_minimal()
#ggsave(native_EU_rel, file="native_EU_rel.png", bg="white")

if(intra_analysis){
  country_native_intra_rel<- country_neophyte_relative[country_neophyte_relative$Neophyte=="native_intra",]
  Europe_nt_in<- left_join(medRegions,country_native_intra_rel,  by= c("Region"="Region"))
  native_intra_EU_rel<-ggplot()+
    geom_sf(data= Europe_nt_in, aes(fill= relative_proportion))+
    scale_fill_viridis_c(option = "magma",begin = 0.1)+
    labs(title = "Native Distribution in European Regions alien elsewhere") +
    theme_minimal()
  #ggsave(native_intra_EU_rel, file="native_intra_EU_rel.png", bg="white")
}


#### 5 CHECK GLONAF ####
###### 5.1 Regions and List #####
# Region shapefile
glonafRegions <- read_sf("../GloNAF_Shapefile", "regions2")
# Make WGS84
glonafRegions <- st_transform(glonafRegions, CRS("+proj=longlat +datum=WGS84"))

# Read list all regions
glonafRegionList<- read.csv("../GloNAF_Shapefile/Region_GloNAF_vanKleunenetal2018Ecology.csv")
# Subset to Europe
glonafRegionList<- glonafRegionList[glonafRegionList$tdwg1_name=="Europe", ]  

# Subset Regions to only those where the ID matches one of the regions in Europe
newGlonaf<- subset(glonafRegions, glonafRegions$OBJIDsic %in% glonafRegionList$OBJIDsic)
glonafRegions<- newGlonaf
# Check whether we have the same number ID in both datasets
all.equal(sort(glonafRegionList$OBJIDsic), sort(glonafRegions$OBJIDsic))
# Merge 
newGlonaf<- left_join(newGlonaf, glonafRegionList, by= "OBJIDsic")




###### 5.2 Species ######
# Read in species
glonafSpecies<-readxl::read_excel('../GloNAF_Shapefile/Taxon_x_List_GloNAF_vanKleunenetal2018Ecology.xlsx')
# Downsize for ease later
glonafSpecies<- glonafSpecies[, c("standardized_name", "region_id")]
# Only Species present in places with region id
glonafSpecies <- glonafSpecies[glonafSpecies$region_id %in% glonafRegionList$region_id, ]
# Check 
all.equal(sort(unique(glonafSpecies$region_id)), sort(glonafRegionList$region_id))

# Combine Species with Region list
glonafSpecies <- left_join(glonafSpecies,glonafRegionList, by="region_id")
glonafSpecies<- glonafSpecies[, c("standardized_name","region_id","tdwg4_name")]
names(glonafSpecies)[names(glonafSpecies) == 'tdwg4_name'] <- 'Region'

# To be able to compare GLONAF and our data, we check the names
setdiff(unique(glonafSpecies$Region), unique(eva_country_neophyte$Region))
setdiff(unique(eva_country_neophyte$Region), unique(glonafSpecies$Region))

# We will take the names of GLONAF for this check-up
correctCountries<- data.frame(Med=c("Rf.NW", "Rf.N","Rf.E", "Rf.C","Rf.S","Rf.K", "Rf.CS","Luxemburg","Bosnia.Herzegovina", "Italy" ,"Czech.Republic", "Greece", "France", 
                                  "Moldavia", "Republic.of.Ireland", "United.Kingdom", "Corsica", "Sicily", "Sardinia", 
                                    "Spain mainland", "Portugal mainland", "Faroes"), 
                              WF=c("Northwest European Russia", "North European Russia", "East European Russia", "Central European Russia", 
                                   "South European Russia", "Kalingrad Region","Russian Caucasia",
                                   "Luxembourg", "Bosnia and Herzegovina", "Italy", "Czech Republic", "Greece",
                                   "France", "Moldova", "Ireland","Great Britain", "Corse", 
                                   "Sicilia", "Sardegna", "Spain", "Portugal", "F?royar"))

# Change our names: first take the index
index <- eva_country_neophyte$Region %in% correctCountries$Med
# Change to correct for each Region
eva_country_neophyte$Region[index] <- correctCountries$WF[match(eva_country_neophyte$Region[index],correctCountries$Med)]
# Quite some countries are not present in GLONAF, we will just check countries we have
setdiff(unique(eva_country_neophyte$Region), unique(glonafSpecies$Region))

# Get names of regions
glonafRegionNames<- unique(glonafSpecies$Region)
MEDRegionNames <- unique(eva_country_neophyte$Region)

# Subset datasets to be able to compare them
glonafSpecies<- glonafSpecies[glonafSpecies$Region %in% MEDRegionNames,]
medSpecies<- eva_country_neophyte[eva_country_neophyte$Region %in% glonafRegionNames,]

# Take only neophytes
medSpecies <- medSpecies[(medSpecies$Neophyte=="extra"|medSpecies$Neophyte=="intra"), ]
# Take only all unique combinations of Region, species and definition.
medUnique<- medSpecies |> group_by(Region, species, Neophyte) |> summarise(n=n())

###### 5.3 Compare #####
# Get all species from GLONAF that are in our database
species_not_MED<- data_frame(country= c(), notMed=c())
species_not_glonaf<- data_frame(country= c(), notGlonaf= c())

# test which are not in dataset
for (i in glonafRegionNames){
  notMed<- setdiff(glonafSpecies$standardized_name[glonafSpecies$Region==i][glonafSpecies$standardized_name[glonafSpecies$Region==i] %in% eva_country_neophyte$species[eva_country_neophyte$Region==i]], medUnique$species[medUnique$Region==i])
  notGlonaf<- setdiff(medUnique$species[medUnique$Region==i],glonafSpecies$standardized_name[glonafSpecies$Region==i][glonafSpecies$standardized_name[glonafSpecies$Region==i] %in% eva_country_neophyte$species[eva_country_neophyte$Region==i]])
  notGlonaf<- cbind(rep(i, times=length(notGlonaf)), notGlonaf)
  notMed<-  cbind(rep(i, times=length(notMed)), notMed)
  species_not_MED<- rbind(species_not_MED,notMed)
  species_not_glonaf<- rbind(species_not_glonaf, notGlonaf)
}

# Quite some species are defined as alien in our database while this is not true in Glonaf (2262 species over 41 countries)
# Maybe worse is that quite some species are not defined as alien in our database but yes in Glonaf (1569 species over 41 countries)
# Randomly checking some species from the latter case against KEW and CABI gives more trust to our own classification


###### 5.4 Origin #####
# Now we will check the origin using the data of Zhang et al 2023 (which is based on the GIFT database)
native<-read.csv("../Neophyte Assignments/wf_dat.csv")
# Remove all _ and make spaces
native$sp_tpl <-gsub("_"," ",native$sp_tpl)

# We make a vector of all unique species in Europe
EUROPE<- unique(native$sp_tpl[native$botanical_continent=="EUROPE"])
NOT_EU<- unique(native$sp_tpl[!native$botanical_continent=="EUROPE"])
  
# We check how much and which species are present in Europe but were defined as from extra European origin
length(which(EUROPE %in% extra_EU))
EUROPE[EUROPE %in% extra_EU]
# These are quite some... I took about 10 random species (5%) and searched their native range, all of which are from extra European origin

# Check how many do come from other continents
length(which(NOT_EU %in% extra_EU))

# Check how many intra european species are found in the dataset
length(which(EUROPE %in% intra_EU))
length(which(NOT_EU %in% intra_EU))
NOT_EU[NOT_EU %in% intra_EU]
# Again quite some, after checking 10 species (5%), only one was not found native to Europe