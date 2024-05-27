
rm(list=ls())

# Load eva and plot data
eva <- read_csv("fullPlotEva_cover_all_layer_cleaned.csv", show_col_types = FALSE)
fullPlotData <- read_csv("fullPlotData_cover_all_layer_cleaned.csv", show_col_types = FALSE)

# Downsample by a factor of 100 if fast is selected
fast <- F
if(fast) {
  fullPlotData <- fullPlotData[runif(length(fullPlotData$PlotObservationID)) > 0.90,]
  eva <- eva[eva$PlotObservationID %in% fullPlotData$PlotObservationID,]
}


#### 1 Cultivated ####

uniqueSpecies<- unique(eva$species)

uniqueSpecies<- as.data.frame(uniqueSpecies)

ecocrop<-readxl::read_excel("C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/Core authors/ECOCROPP.xlsx")
ecocrop<- ecocrop[!is.na(ecocrop$Category),]

ecospecies<- (ecocrop$species)

ecospecies<- sort(uniqueSpecies[uniqueSpecies %in% ecospecies])

test<- eva[eva$species %in% ecospecies,]
test<- c("Zea mays")
test<- eva$PlotObservationID[eva$species %in% test]
test<- test$species[test$`Cover %`>80]
unique(test)
x<-eva[eva$PlotObservationID %in% test,]

unique(test$species)

species_vector <- c("Beta vulgaris", "Agrostis canina", "Agrostis gigantea", "Agrostis stolonifera", 
                    "Agrostis capillaris", "Alopecurus pratensis", "Arrhenatherum elatius", 
                    "Bromus catharticus", "Bromus sitchensis", "Cynodon dactylon", "Dactylis glomerata", 
                    "Festuca arundinacea", "Festuca filiformis", "Festuca ovina", "Festuca pratensis", 
                    "Festuca rubra", "Festuca trachyphylla", "Festulolium", "Lolium multiflorum", 
                    "Lolium perenne", "Lolium x hybridum", "Phalaris aquatica", "Phleum nodosum", 
                    "Phleum pratense", "Poa annua", "Poa nemoralis", "Poa palustris", "Poa pratensis", 
                    "Poa trivialis", "Trisetum flavescens", "Biserrula pelecinus", "Galega orientalis", 
                    "Hedysarum coronarium", "Lathyrus cicera", "Lotus corniculatus", "Lupinus albus", 
                    "Lupinus angustifolius", "Lupinus luteus", "Medicago doliata", "Medicago italica", 
                    "Medicago littoralis", "Medicago lupulina", "Medicago murex", "Medicago polymorpha", 
                    "Medicago rugosa", "Medicago sativa", "Medicago scutellata", "Medicago truncatula", 
                    "Medicago x varia", "Onobrychis viciifolia", "Ornithopus compressus", 
                    "Ornithopus sativus", "Pisum sativum", "Trifolium alexandrinum", "Trifolium fragiferum", 
                    "Trifolium glanduliferum", "Trifolium hirtum", "Trifolium hybridum", 
                    "Trifolium incarnatum", "Trifolium isthmocarpum", "Trifolium michelianum", 
                    "Trifolium pratense", "Trifolium repens", "Trifolium resupinatum", 
                    "Trifolium squarrosum", "Trifolium subterraneum", "Trifolium vesiculosum", 
                    "Trigonella foenum-graecum", "Vicia benghalensis", "Vicia faba", "Vicia pannonica", 
                    "Vicia sativa", "Vicia villosa", "Brassica napus", "Brassica oleracea", 
                    "Phacelia tanacetifolia", "Plantago lanceolata", "Raphanus sativus", 
                    "Arachis hypogea", "Brassica rapa", "Brassica juncea", "Brassica nigra", 
                    "Cannabis sativa", "Carthamus tinctorius", "Carum carvi", "Gossypium spp.", 
                    "Helianthus annuus", "Linum usitatissimum", "Papaver somniferum", "Sinapis alba", 
                    "Glycine max", "Avena nuda", "Avena sativa", "Avena strigosa", "Hordeum vulgare", 
                    "Oryza sativa", "Phalaris canariensis", "Secale cereale", "Sorghum bicolor", 
                    "Triticosecale", "Triticum aestivum", "Triticum turgidum", "Triticum aestivum", 
                    "Zea mays", "Solanum tuberosum")

uniqueSpecies<- unique(eva$species)
agro<- uniqueSpecies[uniqueSpecies %in% species_vector]

# Data on which species are neophytes
native_intra_analysis=F
if(native_intra_analysis){
  species_country_status <- read_csv("eva2_country_status_new.csv")
} else{
  species_country_status<- read_csv("species_country_status_new.csv")
  # or if we want to do it with intra seperately
}

# Assign these names to the eva list
extra_EU <- unique(species_country_status$species[species_country_status$Neophyte=="extra"])
intra_EU <- unique(species_country_status$species[species_country_status$Neophyte=="intra"])
native_intra <- unique(species_country_status$species[species_country_status$Neophyte=="native_intra"])



remove<- read_csv("not_defined.csv",show_col_types = FALSE)
remove<- as.vector(unlist(remove))
eva<- eva[!(eva$species %in% remove),]

# Only observation ID and Region
fullPlot2<- fullPlotData[,c("PlotObservationID","Region")]
# Right join to keep only species present in fullplot (otherwise a lot of NAs)
eva<- right_join(eva, fullPlot2, by = c("PlotObservationID"="PlotObservationID"))
# Join eva and classification
eva_country_neophyte<- left_join(eva, species_country_status, 
                                 by= c("Region"= "Region", "species"= "species"))
eva<- eva_country_neophyte
# Look at how much from every type are present
tt<- as.data.frame(table(eva_country_neophyte$Neophyte))
tt

rm(fullPlot2, eva_country_neophyte, remove)

test <- eva |> group_by(species, Neophyte) |> summarise(n=n())

test <- test[test$species %in% agro, ]

write.csv(test, "cultivated_species_status.csv")

obs<- eva$PlotObservationID[eva$species %in% agro]
test<- fullPlotData[fullPlotData$PlotObservationID %in% obs,]
hist(test$numberOfVascularPlantSpecies[test$numberOfVascularPlantSpecies<=5])
obs_SR<- test$PlotObservationID[test$numberOfVascularPlantSpecies<=5]
test<- eva[eva$PlotObservationID %in% obs_SR & eva$species %in% agro,]
hist(test$`Cover %`)
test<- test[test$`Cover %`>60,]
unique(test$species)


# check overlap two databases
eco<- uniqueSpecies[uniqueSpecies %in% ecospecies]
agro[agro %in% eco]

# check species in agro but  not in the ecocrop database
not_agro<- agro[!agro %in% eco]
length(eva$species[eva$species %in% not_agro])

# check species in the ecocrop database not in the EU agro database
not_eco<- eco[!eco %in% agro]
length(eva$species[eva$species %in% not_eco])


#### 2 Invasives ####

gisd <- read.csv("../EIVE Data/GISD.csv", sep = ";")
gisd_species <- gisd$Species

gisd_species[gisd_species %in% uniqueSpecies]

test <- read.csv("coverClassImpactForCandidates_cover_all_layer_cleaned.csv")

test$invasive <- test$taxa %in% gisd$Species
