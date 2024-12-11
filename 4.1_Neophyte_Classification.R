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
library(readxl)







##### 2 DATA ####
###### 2.1 neophytes ######
# new file created with all alien species in europe (extra+ intra)
# old one from irena
neophyte<- readxl::read_excel("../Neophyte-Impact/country_species-2024-01-30-IA-VK.xlsx", sheet="country_speciesWillem")
# my own one, only species also present in our plots --> all known
neophyte <- read_csv("neophyte_euro.csv", show_col_types = FALSE)
# newest version
neophyte <- read_csv("neophyte_ESy.csv", show_col_types = FALSE)


# compare with Veronika
veronika_aliens <- read_xlsx("willem-alien-data-checked-2024-10-10.xlsx", sheet= "aliens")
veronika_aliens <- veronika_aliens[veronika_aliens$Status != veronika_aliens$Status.checked.VK,]
veronika_aliens <- veronika_aliens[!is.na(veronika_aliens$Rationale),]


# glonaf
veronika_glonaf <- read_xlsx("willem-alien-data-checked-2024-10-10.xlsx", sheet= "species-not-GloNAF")
veronika_med <- read_xlsx("willem-alien-data-checked-2024-10-10.xlsx", sheet= "species-not-MED")
veronika_med <- veronika_med[!veronika_med$reasoning=="POWO,Euro+Med = native AND/OR checklist = archaeophyte, uncertain",]

# how many species in the newest version
test <- veronika_aliens[veronika_aliens$`Aggregated name` %in% veronika_med$notMed,]


# problems
veronika_problems <- read_xlsx("willem-alien-data-checked-2024-10-10.xlsx", sheet= "merging-problems-status")
veronika_problems <- veronika_problems[!veronika_problems$Conflict=="no",]

###### 2.2 Eva and Header #####
# Load and make smaller to conocate
fullPlotEva <- read_csv("fullPlotEva_ESy.csv", show_col_types = FALSE)
eva2<- fullPlotEva[,c("PlotObservationID","species", "irena","Matched concept", "name")]
fullPlotData<- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
fullPlot2<- fullPlotData[,c("PlotObservationID","Region")]


# Load summarized data
# these is the full region-species database
# these include unknown classifications
species_country <- read_csv("country_species_ESy.csv", show_col_types = FALSE)

# Join eva and header (species per plot and plot info respectively)
eva_country<- left_join(eva2, fullPlot2)
# Join eva_country with neophyte data from Irena --> gives all alien species per country
eva_country_neophyte<- left_join(eva_country[,-c(2:4)], species_country, by= c("Region"= "Region", "name"="name"))
eva_country_neophyte<- subset(eva_country_neophyte, select=-c(n, SeqID))
# check columns
head(eva_country_neophyte)
colnames(eva_country_neophyte)[8]<- "Neophyte"
# species = our new species name
# irena = name irena file
# matched concept = original name
# name = name ESy 



# Take dataset with all neophytes within EU
neophyteDefEU<- (eva_country_neophyte[(!(eva_country_neophyte$Neophyte=="native")& !is.na(eva_country_neophyte$Neophyte)),])
any(is.na(eva_country_neophyte$Neophyte ) & eva_country_neophyte$name != 'Plant')
any(is.na(eva_country_neophyte$Neophyte ) & eva_country_neophyte$name == 'Plant')


# Retrieve all names of all alien species to Europe
neophyteNamesEU<- unique(neophyteDefEU$species)
#In total, we have 1772 alien species in European countries (before removing NA dates, after doing so 1686, (1699)) -->1661 after ESy


# Check whether joining worked fine
Poland <- unique(sort(neophyteDefEU$name[neophyteDefEU$Region=="Poland"]))
Poland2 <-(unique(neophyte$name[neophyte$Region== "Poland"& neophyte$statusNew!="native"]))
setdiff(Poland, Poland2)
setdiff(Poland2, Poland)


colnames(species_country)[9] <- "Neophyte"

# change some species with Veronika
for(i in 1: nrow(veronika_aliens)){
  species_country$Neophyte[species_country$Region==veronika_aliens$Region[i] & species_country$name==veronika_aliens$`Accepted name`[i]] <- veronika_aliens$Status.checked.VK[i]
}
i=1
# change some species with Veronika
for(i in 1: nrow(veronika_problems)){
  species_country$Neophyte[species_country$Region==veronika_problems$Region[i] & species_country$name==veronika_problems$name[i]] <- veronika_problems$`Aggregated.status.VK (with respect to a given region)`[i]
}


#write_csv(species_country, "country_species_ESy.csv")


intra_analysis=T
if(intra_analysis){
  # Check whether some intra_EU species are native in other regions
  native_names<-unique(species_country$species[species_country$Neophyte=="native"])
  intra_EU <- unique(species_country$species[species_country$Neophyte=="intra"])
  extra_EU <- unique(species_country$species[species_country$Neophyte=="extra"])
  
  sum(native_names %in% intra_EU) #818 (before around 848) (823) --> 789
  native_intra<- native_names[native_names %in% intra_EU]
  
  
  # Here we make a new dataframe with those species that are intra in a region in europe as native_intra
  # to check whether the effect is just species dependent
  eva2<- eva_country_neophyte
  eva2$Neophyte[eva2$Neophyte=="native" & eva2$species %in% native_intra]<-"native_intra"
  length(unique(eva2$species[eva2$Neophyte=="extra"]))
  length(unique(eva2$species[eva2$Neophyte=="intra"]))
  length(unique(eva2$species[eva2$Neophyte=="native"]))  
  length(unique(eva2$species[eva2$Neophyte=="native_intra"])) 
  
  # check whether the classification was successful for all species
  sum(unique(eva2$species[eva2$Neophyte=="native"]) %in% intra_EU)
  
  country_neophyte <- eva2 |> distinct(Region, Neophyte, species) %>% group_by(Region, Neophyte) |> summarise(n=n())
  table(country_neophyte['Neophyte'])
  eva_country_neophyte<- eva2
  
  species_country$Neophyte[species_country$Neophyte=="native" & species_country$species %in% native_intra]<-
    "native_intra"
  species_country_status<- read_csv("country_species_ESy.csv", show_col_types = FALSE)
}

###### 2.3 Crops ######
species_data <- data.frame(old = c( "crop vineyard", "crop barley", "crop maize", "crop wheat", "crop potato","crop hops", "crop pea", "crop digitalis", "crop tomato", 
                                    "crop currant","crop cabbage", "crop bean", "crop onion", "crop asparagus", "crop spelt","crop carrot", "crop rhubarb", "crop clover", 
                                    "crop tobacco", "crop zucchini", "crop parsley", "crop gladiolus", "crop salsify", "crop radish", "crop lettuce","crop buckwheat"),
                           new = c("Vitis vinifera", "Hordeum vulgare", "Zea mays", "Triticum aestivum", "Solanum tuberosum", "Humulus lupulus", "Pisum sativum", 
                                   "Digitalis purpurea", "Solanum lycopersicum", "Ribes rubrum","Brassica oleracea", "Phaseolus vulgaris", "Allium cepa", 
                                   "Asparagus officinalis", "Triticum aestivum","Daucus carota", "Rheum rhabarbarum", "Trifolium repens", "Nicotiana tabacum", 
                                   "Cucurbita pepo", "Petroselinum crispum", "Gladiolus grandiflorus", "Tragopogon porrifolius", "Raphanus sativus", "Lactuca sativa",
                                   "Fagopyrum esculentum"))


species_country$name[species_country$name %in% species_data$old] <- species_data$new[match(species_country$name[species_country$name %in% species_data$old], species_data$old)]

# check for duplicats
dup <- species_country[duplicated(species_country[,c(1:5)]) |duplicated(species_country[,c(1:5)], fromLast=T), ]

#write_csv(species_country, "country_species_ESy.csv")

# the remainder is not necessary anymore, we have delineated all species
further= F
if(further){

species_country_status<- read_csv("country_species_ESy.csv", show_col_types = FALSE)
aliens <- species_country_status[species_country_status$Neophyte=="intra"| species_country_status$Neophyte=="extra",c(1:2,4:5, 9) ] 
colnames(aliens)<- c("Region","Accepted name","Original EVA name","Aggregated name","Status")
write_csv(aliens, "aliens.csv")  
  
###### 2.4 extra-EU #######
# Data on which species are neophytes from outside of Europe
neophyteDefinitions <- read_csv("../Neophyte-Impact/Neophyte Assignments/UniqueTaxaEurope-2023-04-23.csv", show_col_types = FALSE)
# get names eva
eva_names <- unique(eva2[, c(2:4)])


# change names to all accepted
# first match original concept with the accepted species name
neophyteDefinitions$name <- eva_names$species[match(neophyteDefinitions$Matched.concept, eva_names$species)]
# 6305 species not yet--> match with irena
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$Matched.concept[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$irena)]
# 4686 species not yet--> match matched concept
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$Matched.concept[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$`Matched concept`)]
# 2793 species not yet--> link the name in neophyteDefinitions (hence given by irena)
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$species[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$species)]
# 2372 not yet
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$species[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$irena)]
# 2366 species not yet
neophyteDefinitions$name[is.na(neophyteDefinitions$name)] <- eva_names$species[match(neophyteDefinitions$species[is.na(neophyteDefinitions$name)], 
                                                                                     eva_names$`Matched concept`)]
# check species without a match, likely removed from our data so no problem
no_name<- neophyteDefinitions[is.na(neophyteDefinitions$name),]


# redefine
neophyteDefinitions <- neophyteDefinitions[, c(7, 4,6,3)]
# species for our accepted name, 'name' for the original name in the dataset
colnames(neophyteDefinitions)<- c("species","name","neophyte","exclude")
# check duplicates
neophyteDefinitions <- neophyteDefinitions[!duplicated(neophyteDefinitions),]


# we check again whether our change was successful
# check for all duplicates in our accepted species names (from last true as well to also obtain the other duplicates)
x <- neophyteDefinitions[(duplicated(neophyteDefinitions[,c(1)])| duplicated(neophyteDefinitions[,c(1)], fromLast=TRUE)),]
# check whether all neophyte classifications overlap (this is the most crucial)
x <- x[!(duplicated(x[,c(1,3)]) | duplicated(x[,c(1,3)], fromLast=TRUE)),]


# has to be changed for some species
change<- data.frame(species= c("Avena sativa","Cynara scolymus", "Fragaria moschata","Gnaphalium",
                               "Raphanus raphanistrum subsp. raphanistrum","Rubus canadensis"),
                    status= c("arch","native","native","native","native","neo"))


# assign to neophyte
neophyteDefinitions$neophyte[neophyteDefinitions$species %in% change$species] <- 
  change$status[match(neophyteDefinitions$species[neophyteDefinitions$species %in% change$species],change$species)]
# we check again whether our change was successful
x <- neophyteDefinitions[(duplicated(neophyteDefinitions[,c(1)])| duplicated(neophyteDefinitions[,c(1)], fromLast=TRUE)),]
x <- x[!(duplicated(x[,c(1,3)]) | duplicated(x[,c(1,3)], fromLast=TRUE)),]
# yes, now all are true


# remove all complete duplicates
# this can be done safely as we know that all duplicates in neophytedefinitions (which does not contain regions) do now contain the same info
neophyteDefinitions <- neophyteDefinitions[!(duplicated(neophyteDefinitions)),]


# check duplicates of species definitions --> some species might be duplicated but the irena or matched concept not
dup<- neophyteDefinitions[duplicated(neophyteDefinitions$species)| duplicated(neophyteDefinitions$species, fromLast=TRUE),]
# there are 2330 NA values in this file --> weird as there are 2366 species with no name
sum(is.na(dup$species))
# check whether all these species are in no_name (this file was created before our name change so the old name is still species)
nrow(dup[dup$name %in% no_name$species,])
# all species here are also in the file
# they are also almost all in neophyteDefinitions
nrow(no_name[no_name$species %in% dup$name,])


# there are also some non-NA values in the file, check
dup[!is.na(dup$species),]
# these are species differing between irena and us
# check whether for them all alien classifications applies
x <- dup[!(duplicated(dup[,c(1,3)]) | duplicated(dup[,c(1,3)], fromLast=TRUE)),]
# these can safely be removed


# remove
neophyteDefinitions <- neophyteDefinitions[!duplicated(neophyteDefinitions$species),]
# remove the one species still na
neophyteDefinitions<- neophyteDefinitions[!is.na(neophyteDefinitions$species),]


# join eva names (all species present) and neophyte definitions (all species alien to Europe, cleaned so that it only contains species present)
species <- left_join(eva_names, neophyteDefinitions, by= c("species"="species"))
# check whether it helps to also do this for Irena
species[is.na(species$neophyte),] <- left_join(species[is.na(species$neophyte),-c(4:6)], neophyteDefinitions, by= c("irena"="species"))
# 196 --> same as previously
unique(species$species[is.na(species$neophyte)])


# Some species should be excluded --> Check which 
exclude <- neophyteDefinitions$species[neophyteDefinitions$neophyte == "exclude"]
exclude <- exclude[exclude %in% fullPlotEva$species]
# Maybe best to remove these from eva --> otherwise bias in application of names


# But first check how they are incorporated in the new file --> not neophyte
# remove from intra European alien file
neophyteDefEU <- neophyteDefEU[!(neophyteDefEU$species %in% exclude),]
neophyteNamesEU <- neophyteNamesEU[!(neophyteNamesEU %in% exclude)]
# all eva observations --> two observations removed
eva_country_neophyte <- eva_country_neophyte[!(eva_country_neophyte$species %in% exclude),]
# remove from extra european alien file --> removed one
neophyteDefinitions <- neophyteDefinitions[!(neophyteDefinitions$species %in% exclude),]
neophyteNames <- neophyteDefinitions$species[neophyteDefinitions$neophyte=="neo"]


# Some species are archeophytes
arch <- neophyteDefinitions$species[neophyteDefinitions$neophyte == "arch"]
arch <- arch[arch %in% fullPlotEva$species]
# We will not do anything with them now, but it is possible to adjust the code in the end to also check archeophytes 
# and whether they different (85, 70)


###### 2.5 SUMMARY #####
# All extra EU neophytes defined neo
# dataframe (neophyteDefinitions)
neophyteNames
# All aliens and the region they are present
# dataframe (neophyte)
neophyteNamesEU



#### 3 ANALYSIS ####
###### 3.1 intra-EU ######
# Take names all intra-EU species --> also possible (setdiff(neophyteNamesEU, neophyteNames)) (899)
# these are all species that were alien in euroe (intra + extra) minus the extra european alien
intra_EU<- neophyteNamesEU[!(neophyteNamesEU %in% neophyteNames)]


# The sums do not add up to all species in neophyteNamesEU, we check which species are present in neophyteNames and not in neophyteNamesEU (1723)
all<- c(intra_EU, neophyteNames)
not_defined<- setdiff(all, neophyteNamesEU) # we checked other way around, is empty character vector: setdiff(neophyteNamesEU, all)
# CONCLUSION: there are 25 (21, 24) species that were defined extra-European in the old dataset but were now included in the list as native, causing the
# mismatch between both datasets. 


# Extra test to see where they can be found
country_species_number<- eva_country_neophyte |> group_by(Region, species, Neophyte) |> summarise( n= n())
print(country_species_number[country_species_number$species %in% not_defined, ], n= 45)
# most species are found in Turkey --> in old one defined as extra-European but now included.
test<- (neophyte[neophyte$species %in% not_defined,])
test<- test[!(test$species=="Cephalophysis species"),]
# write.csv(test, 'species_not_defined.csv', row.names=FALSE)


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


# add these species to extra
add_extra<- as.character(c("Citrus x limon","Citrus × limon","Cucumis melo","Phoenix dactylifera"))
# add species to intra
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
# eva_country_neophyte$Neophyte[is.na(eva_country_neophyte$Neophyte)] <- "native"
eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% extra_EU] <- "extra"
eva_country_neophyte$Neophyte[eva_country_neophyte$Neophyte=="neo"]<- "intra"

# unknowns 257 species  
length(unique(eva_country_neophyte$species[is.na(eva_country_neophyte$Neophyte)]))
# total of 3434 observations
unknown <- (eva_country_neophyte[is.na(eva_country_neophyte$Neophyte),])
# total of 269 species country relations
unknown <- unknown[!duplicated(unknown[, c(2:8)]),]
# some of these are in neophyte but not in our database
unknown$species[!unknown$species %in% neophyte$species]


# Check archeophytes
unique(eva_country_neophyte$Neophyte[eva_country_neophyte$species %in% arch])
# one is intra-european
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
# summarise again based on region and species --> this time we have extra and intra instead of neo
country_status <- eva_country_neophyte |> group_by(Region, species, irena, `Matched concept`, Neophyte) |> summarise(n=n())


# check difference with previously made data --> we removed one species (see exclude)
setdiff(country_species_number[, c(1:2)], country_status[, c(1:2)] )


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
  length(unique(eva2$species[eva2$Neophyte=="native"])) # 18794
  length(unique(eva2$species[eva2$Neophyte=="native_intra"])) # 824 (how come one more?)
  
  # check whether the classification was successful for all species
  sum(unique(eva2$species[eva2$Neophyte=="native"]) %in% intra_EU)
  
  country_neophyte <- eva2 |> distinct(Region, Neophyte, species) %>% group_by(Region, Neophyte) |> summarise(n=n())
  table(country_neophyte['Neophyte'])
  eva_country_neophyte<- eva2
}




###### 3.3 native SR #####
# aggregate and count number of species
aggregatedEVA <- eva_country_neophyte |>  group_by(PlotObservationID, Neophyte) |>  summarise(numberOfVascularPlantSpecies = n())


# Only native SR is relevant here (does not make sense to look at the influence of alien species on alien species)
nativeSpR <- aggregatedEVA[(aggregatedEVA$Neophyte=="native" | aggregatedEVA$Neophyte=="native_intra") & !is.na(aggregatedEVA$Neophyte),]
nativeSpR <- subset(nativeSpR, select= -c(Neophyte))
names(nativeSpR)[names(nativeSpR)=="numberOfVascularPlantSpecies"]<- "nativeSR"
# add all per plot (as native and native_intra both are native species)
nativeSpR <- nativeSpR |> group_by(PlotObservationID) |> summarise(nativeSR = sum(nativeSR))

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
  # there is a difference with extra and intra --> 3437 difference
  table(eva_country_neophyte['Neophyte'])
  sum(eva_country_neophyte$Neophyte=="extra", na.rm=T)+sum(eva_country_neophyte$Neophyte=="intra", na.rm=T)
}
 
# check with natives as well
sum(fullPlotData$numberOfVascularPlantSpecies, na.rm=T)
# check sum per species
sum(eva_country_neophyte$Neophyte=="native", na.rm=T)+ 
  sum(eva_country_neophyte$Neophyte=="extra"| eva_country_neophyte$Neophyte=="intra" | eva_country_neophyte$Neophyte=="native_intra", na.rm=T)
# check currently na
sum(is.na(eva_country_neophyte$Neophyte))
# the difference of 3 persists  


# difference is also present in difference eva_country and eva_country_neophyte --> check!
# difference is caused by the removal of species with the exclude and not defined part --> check where it is present and reduce with 1
remove<- c(not_defined, exclude)
if(nativeSR_calculation){
  # species present in 3 sites
  remove_observations <- eva_country$PlotObservationID[eva_country$species %in% remove]
  # reduce species richness here with 1
  fullPlotData$numberOfVascularPlantSpecies[fullPlotData$PlotObservationID %in% remove_observations]<-
    fullPlotData$numberOfVascularPlantSpecies[fullPlotData$PlotObservationID %in% remove_observations]-1
}
# checked --> correct



###### 3.4 Save #####
# save the species x country x neophyte status file
if(intra_analysis){
  eva2_country_status<- eva2 |> group_by(Region, species, Neophyte) |> summarise(n=n())
  eva2_country_status<- eva2_country_status[,-4]
  #write.csv(eva2_country_status,"eva2_country_status_euro.csv", row.names = FALSE)
  } else {
    species_country_status<- eva_country_neophyte |> group_by(Region, species, Neophyte) |> summarise(n=n())
    species_country_status<- species_country_status[,-4]
    #write.csv(species_country_status,"species_country_status_new.csv", row.names = FALSE)
}

# save the change in the fullPlotData file (removed 3 species in number of vascular plant species)
#write_csv(fullPlotData, "fullPlotData_euro.csv")

# remove species
remove<- c(not_defined, exclude)
#write.csv(remove, "not_defined.csv", row.names=FALSE)

# made last time with all species in euro
# write_csv(new_names, "new_names.csv")
#write_csv(country_status, "species_status_region.csv")

# change species check for Koenraad
#country_status$Neophyte[country_status$species=="Origanum vulgare" & country_status$Region=="Belgium"] <- 'native'
}

###### 3.5 Koenraad analysis #####
future<- read_excel("~/Boeren/Impact/eva_neophytes/Core authors/Species_FutureNature.xlsx")
future$Region <- c("Belgium")
country_status <- read_csv("species_status_region.csv",show_col_types = FALSE)

future <- left_join(future, country_status, by= c("Species"="species", "Region"="Region"))

future[is.na(future$Neophyte), ] <- left_join(future[is.na(future$Neophyte), c(1:2)], 
                                              country_status, by= c("Species"="irena", "Region"="Region"))

future$Method<-NA
future$Method[!is.na(future$Neophyte)] <- "Irena"


neophyte <- read_csv("neophyte_euro.csv", show_col_types = FALSE)
neophyteDefinitions <- read_csv("../Neophyte-Impact/Neophyte Assignments/UniqueTaxaEurope-2023-04-23.csv", show_col_types = FALSE)

missing <- future[is.na(future$Neophyte),]


missing$Neophyte[missing$Species %in% neophyteDefinitions$species] <- neophyteDefinitions$statusEurope[neophyteDefinitions$species %in% missing$Species]

missing <- left_join(missing, neophyteDefinitions[, c(2,4,6)], by= c("Species"="species"))

missing$Species[!missing$Species %in% neophyteDefinitions$species]


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

glonafSpecies <- glonafSpecies[glonafSpecies$Region=="Belgium",]




missing <- left_join(missing, glonafSpecies, by= c("Species"= "standardized_name"))

future$Genus <- gsub("([A-Za-z]+).*", "\\1", future$Species)



library(kewr)
# create list to add our results to
test <- list()
#create vector to store things deemed incorrect (so species in our database alien and not in glonaf)


y <- c()
z<- c()
for( i in 1: nrow(future)){
  # search
  x<- search_powo(list(name=future[i,1], distribution= future[i,2]))
  x$results
  tmp <- tidy(x)
  # add to list
  test[[i]] <- x
  # if the result returns a species this means that it is recognized and present in the country so we should check our classification
  if(nrow(tmp)>0){
    if(ncol(tmp)<10){
      tmp$images<- NA
    }
    # add name of the species
    tmp <- cbind(tmp,future[i,1])
    # store
    y <- rbind(y, tmp)
    z <- rbind(z, y[y$name %in% future[i,1],])
  }
}


full <- future[future$Species%in% y$Species,]
full$Method[is.na(full$Method)] <- "POWO"

intra <- c("Lathyrus aphaca","Jacobaea maritima","Galatella sedifolia",
           "Digitaria sanguinalis","Briza minor","Brachypodium phoenicoides")

full$Neophyte[full$Species %in% intra] <- "intra"
full$Neophyte[is.na(full$Neophyte)] <- "native"

not_full <-future[!future$Species%in% y$Species,]

not_full$Neophyte[not_full$Species %in% glonafSpecies$standardized_name] <- "intra"
not_full$Method[not_full$Species %in% glonafSpecies$standardized_name] <- "GLONAF"

length(not_full$Species[not_full$Species %in% neophyteDefinitions$species])
not_full[not_full$Species %in% neophyteDefinitions$species[neophyteDefinitions$statusEurope=="neo"],]

not_full$Neophyte[not_full$Species %in% neophyteDefinitions$species] <- 'intra'
not_full$Method[not_full$Species %in% neophyteDefinitions$species] <- 'not extra'


not_full$Neophyte[!not_full$Species %in% neophyteDefinitions$species] <-  c("intra","intra","native")  
not_full$Method[!not_full$Species %in% neophyteDefinitions$species] <-  "own"


Species_FutureNature <- rbind(full, not_full)
Species_FutureNature <- Species_FutureNature[, c(1,2,5,7)]


##### 4 MAP ####
###### 4.1 MED regions #####
# We load the Med regions and give correct CRS
medRegions <- read_sf("../Europe-regions-shapefiles-2023", "Emed_regions")
medRegions <- st_transform(medRegions, CRS("+proj=longlat +datum=WGS84"))

# join with the number of species per class per area
density <- left_join(medRegions, country_neophyte, by= c("Region"="Region"))
# calculate the density
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
extra_EU_plot
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
#glonafSpecies <- unique(glonafSpecies[, c(1,3:4)])

# To be able to compare GLONAF and our data, we check the names
setdiff(unique(glonafSpecies$Region), unique(species_country_status$Region))
setdiff(unique(species_country_status$Region), unique(glonafSpecies$Region))

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
index <- species_country_status$Region %in% correctCountries$Med
# Change to correct for each Region
species_country_status$Region[index] <- correctCountries$WF[match(species_country_status$Region[index],correctCountries$Med)]
# Quite some countries are not present in GLONAF, we will just check countries we have
setdiff(unique(species_country_status$Region), unique(glonafSpecies$Region))

# Get names of regions
glonafRegionNames<- unique(glonafSpecies$Region)
MEDRegionNames <- unique(species_country_status$Region)

# Subset datasets to be able to compare them
glonafSpecies<- glonafSpecies[glonafSpecies$Region %in% MEDRegionNames,]
medSpecies<- species_country_status[species_country_status$Region %in% glonafRegionNames,]

# Take only neophytes
medSpecies <- medSpecies[(medSpecies$Neophyte=="extra"|medSpecies$Neophyte=="intra"), ]
# Take only all unique combinations of Region, species and definition.
medUnique<- medSpecies |> group_by(Region,name, species, Neophyte) |> summarise(n=n())

# give GLONAF our definition of species
# Eva
eva <- read_csv("fullPlotEva_ESy.csv", show_col_types = FALSE)
eva_names <- unique(eva[, c(2:7)])

# link name
glonafSpecies$name <- eva_names$name[match(glonafSpecies$standardized_name, eva_names$name)]
glonafSpecies$name[is.na(glonafSpecies$name)] <- eva_names$name[match(glonafSpecies$standardized_name[is.na(glonafSpecies$name)], eva_names$species)]
glonafSpecies$name[is.na(glonafSpecies$name)] <- eva_names$name[match(glonafSpecies$standardized_name[is.na(glonafSpecies$name)], eva_names$irena)]
glonafSpecies$name[is.na(glonafSpecies$name)] <- eva_names$name[match(glonafSpecies$standardized_name[is.na(glonafSpecies$name)], eva_names$`Matched concept`)]
glonafSpecies$name[is.na(glonafSpecies$name)] <- eva_names$name[match(glonafSpecies$standardized_name[is.na(glonafSpecies$name)], eva_names$`Turboveg2 concept`)]

glonafSpecies <- glonafSpecies[!is.na(glonafSpecies$name),]


CountRegions <- glonafRegionList |> group_by(tdwg4_name) |> summarise(n=n())
glonafSpecies <- glonafSpecies |> group_by(name, standardized_name, Region) |> summarise(n=n())
glonafSpecies$total <- CountRegions$n[match(glonafSpecies$Region, CountRegions$tdwg4_name)]
glonafSpecies$rel <- glonafSpecies$n/glonafSpecies$total

part<- glonafSpecies[(glonafSpecies$rel<1),]
#glonafSpecies <- glonafSpecies[!(glonafSpecies$rel<1),]

###### 5.3 Compare #####
# Get all species from GLONAF that are in our database
species_not_MED<- data_frame(country= c(), notMed=c())
species_not_glonaf<- data_frame(country= c(), notGlonaf= c())

# test which are not in dataset
for (j in 1:length(glonafRegionNames)){
  i <- glonafRegionNames[j]
  notMed<- setdiff(glonafSpecies$name[glonafSpecies$Region==i][glonafSpecies$name[glonafSpecies$Region==i] %in% species_country_status$name[species_country_status$Region==i]], medUnique$name[medUnique$Region==i])
  notGlonaf<- setdiff(medUnique$name[medUnique$Region==i],glonafSpecies$name[glonafSpecies$Region==i][glonafSpecies$name[glonafSpecies$Region==i] %in% species_country_status$name[species_country_status$Region==i]])
  notGlonaf<- cbind(rep(i, times=length(notGlonaf)), notGlonaf)
  notMed<-  cbind(rep(i, times=length(notMed)), notMed)
  # species in glonaf and not in ours
  species_not_MED<- rbind(species_not_MED,notMed)
  # species in ours not in glonaf
  species_not_glonaf<- rbind(species_not_glonaf, notGlonaf)
}

sum(species_not_glonaf$notGlonaf %in% part$name)
sum(species_not_MED$notMed %in% part$name)

write_csv(species_not_MED,"species_not_MED.csv")
write_csv(species_not_glonaf,"species_not_GLONAF.csv")
write_csv(eva_names, "eva_names.csv")

species_MED<- data_frame(country= c(), notMed=c())
species_glonaf<- data_frame(country= c(), notGlonaf= c())

for (j in 1:length(glonafRegionNames)){
  i <- glonafRegionNames[j]
  notMed<- intersect(glonafSpecies$name[glonafSpecies$Region==i][glonafSpecies$name[glonafSpecies$Region==i] %in% species_country_status$name[species_country_status$Region==i]], medUnique$name[medUnique$Region==i])
  notGlonaf<- intersect(medUnique$name[medUnique$Region==i],glonafSpecies$name[glonafSpecies$Region==i][glonafSpecies$name[glonafSpecies$Region==i] %in% species_country_status$name[species_country_status$Region==i]])
  notGlonaf<- cbind(rep(i, times=length(notGlonaf)), notGlonaf)
  notMed<-  cbind(rep(i, times=length(notMed)), notMed)
  # species in glonaf and not in ours
  species_MED<- rbind(species_MED,notMed)
  # species in ours not in glonaf
  species_glonaf<- rbind(species_glonaf, notGlonaf)
}


# Quite some species are defined as alien in our database while this is not true in Glonaf (2331 species over 41 countries)
# Maybe worse is that quite some species are not defined as alien in our database but yes in Glonaf (1568 species over 41 countries)
# Randomly checking some species from the latter case against KEW and CABI gives more trust to our own classification

check <- eva[eva$name=="Origanum vulgare",]
fullPlotData<- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)
check2 <- fullPlotData[fullPlotData$PlotObservationID %in% check$PlotObservationID,]


# check in kew whether the species is present in that region
library(kewr)

# create list to add our results to
test <- list()
#create vector to store things deemed incorrect (so species in our database alien and not in glonaf)
y <- c()
for( i in 1: nrow(species_not_glonaf)){
  # search
  x<- search_powo(list(species=species_not_glonaf[i,2], distribution= species_not_glonaf[i,1]), filters=c("accepted"))
  # create dataframe
  tmp <- tidy(x)
  # add to list
  test[[i]] <- tmp
  # if the result returns a species this means that it is recognized and present in the country so we should check our classification
  if(nrow(tmp)>0){
    if(ncol(tmp)<10){
      tmp$images<- NA
    }
    # add name of the species
    tmp <- cbind(tmp,species_not_glonaf[i,2])
    # store
    y <- rbind(y, tmp)
  }
}
# this seems to be working, although all returned species are incorrect
# all others seem not be present in these regions naturally




# I want to check this also for species present naturally
# read the file with species x region classification
species <- read_csv("species_status_region.csv", show_col_types = FALSE)


# create list to add our results to
test <- list()
#create vector to store things deemed incorrect (so species in our database alien and not in glonaf)
y <- c()
for( i in 1: nrow(species)){
  # search
  x<- search_powo(list(species=species[i,2], distribution= species[i,1]), filters=c("accepted"))
  # create dataframe
  tmp <- tidy(x)
  # add to list
  test[[i]] <- tmp
  # if the result returns a species this means that it is recognized and present in the country so we should check our classification
  if(nrow(tmp)>0){
    if(ncol(tmp)<10){
      tmp$images<- NA
    }
    # add name of the species
    tmp <- cbind(tmp,species[i,2])
    # store
    y <- rbind(y, tmp)
  }
}



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


##### 6 CHECK NSR ####
###### 6.1 Data #####
eva <- read.csv( "country_species_ESy.csv")

# not possible due to lack data here, see: https://bien.nceas.ucsb.edu/bien/tools/nsr/nsr-api/