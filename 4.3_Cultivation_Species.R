
rm(list=ls())
library(readr)
library(dplyr)

# Load eva and plot data
eva <- read_csv("fullPlotEva_ESy_1980.csv", show_col_types = FALSE)
fullPlotData <- read_csv("fullPlotData_ESy_1980.csv", show_col_types = FALSE)

# Downsample by a factor of 100 if fast is selected
fast <- F
if(fast) {
  fullPlotData <- fullPlotData[runif(length(fullPlotData$PlotObservationID)) > 0.90,]
  eva <- eva[eva$PlotObservationID %in% fullPlotData$PlotObservationID,]
}


#### 1 Cultivated ####
# species list to check ecocrop
uniqueSpecies<- unique(eva$name)

uniqueSpecies<- as.data.frame(uniqueSpecies)

ecocrop<-readxl::read_excel("C:/Users/u0166342/Documents/Boeren/Impact/eva_neophytes/Core authors/ECOCROPP.xlsx")
ecocrop<- ecocrop[!is.na(ecocrop$Category),]

ecospecies<- (ecocrop$species)

# this has to be changed if further
eva$name[eva$name=="Medicago sativa"] <- "Medicago sativa aggr."


# change crop names
# change crops 
crops <- eva[grepl("crop ", eva$name, ignore.case = TRUE),]
# Create the dataframe
species_data <- data.frame(old = c( "crop vineyard", "crop barley", "crop maize", "crop wheat", "crop potato","crop hops", "crop pea", "crop digitalis", "crop tomato", 
                                    "crop currant","crop cabbage", "crop bean", "crop onion", "crop asparagus", "crop spelt","crop carrot", "crop rhubarb", "crop clover", 
                                    "crop tobacco", "crop zucchini", "crop parsley", "crop gladiolus", "crop salsify", "crop radish", "crop lettuce","crop buckwheat"),
                          new = c("Vitis vinifera", "Hordeum vulgare", "Zea mays", "Triticum aestivum", "Solanum tuberosum", "Humulus lupulus", "Pisum sativum", 
                                  "Digitalis purpurea", "Solanum lycopersicum", "Ribes rubrum","Brassica oleracea", "Phaseolus vulgaris", "Allium cepa", 
                                  "Asparagus officinalis", "Triticum aestivum","Daucus carota", "Rheum rhabarbarum", "Trifolium repens", "Nicotiana tabacum", 
                                  "Cucurbita pepo", "Petroselinum crispum", "Gladiolus grandiflorus", "Tragopogon porrifolius", "Raphanus sativus", "Lactuca sativa",
                                  "Fagopyrum esculentum"))


eva$name[eva$name %in% species_data$old] <- species_data$new[match(eva$name[eva$name %in% species_data$old], species_data$old)]

# create our own vector
species_vector <- c("Beta vulgaris subsp. vulgaris","Beta vulgaris", "Secale cereale","Triticum aestivum","Pisum sativum", "Solanum tuberosum",
                    "Hordeum vulgare", "Avena sativa", "Zea mays", "Brassica napus","Oryza sativa","Triticum turgidum","Glycine max",
                    "Brassica oleracea", "Sorghum bicolor","Helianthus annuus","Arachis hypogaea","Raphanus sativus", "Medicago sativa aggr.",
                    "Trifolium incarnatum","Phacelia tanacetifolia","Fagopyrum esculentum", "Vitis vinifera", "Hordeum vulgare", "Zea mays", "Triticum aestivum", 
                    "Solanum tuberosum", "Humulus lupulus", "Pisum sativum", 
                    "Digitalis purpurea", "Solanum lycopersicum", "Ribes rubrum","Brassica oleracea", "Phaseolus vulgaris", "Allium cepa", 
                    "Asparagus officinalis", "Triticum aestivum","Daucus carota", "Rheum rhabarbarum", "Trifolium repens", "Nicotiana tabacum", 
                    "Cucurbita pepo", "Petroselinum crispum", "Gladiolus grandiflorus", "Tragopogon porrifolius", "Raphanus sativus", "Lactuca sativa",
                    "Fagopyrum esculentum")
length(unique(species_vector$species))

# check number of observations higher than 50
x <- eva[eva$name %in% species_vector & eva$`Cover %`>=50,]
y<- x |> group_by(name) |> summarise(n=n())

# remove plots from both header and eva
plots_cultivated <- eva$PlotObservationID[eva$name %in% species_vector & eva$`Cover %`>=50]
fullPlotData<- fullPlotData[!fullPlotData$PlotObservationID %in% plots_cultivated, ]                    
eva <- eva[!eva$PlotObservationID %in% plots_cultivated,]

# check status species
species_country_status<- read_csv("country_species_ESy_1980.csv", show_col_types = FALSE)
species_vector <- as.data.frame(species_vector)
colnames(species_vector)<- "species"
species_vector$status <- species_country_status$Neophyte[match(species_vector$species, species_country_status$name)]


#### 2 Changes ####
remove <- c("Entodon concinnus","Hygrohypnum luridum")

which(eva$name %in% remove)

eva <- eva[!eva$name %in% remove,]


#write_csv(eva, "fullPlotEva_ESy_1980.csv")
#write_csv(fullPlotData, "fullPlotData_ESy_1980.csv")
#write_csv(y, "removed_cultivated.csv")


#### 3 Time ####
# here again check time
hist(fullPlotData$Date, breaks= 20)
summary(fullPlotData$Date)

sum(fullPlotData$Date < as.Date("1970-01-01"))/ nrow(fullPlotData)

#### 4 Invasives ####
gisd <- read.csv("../EIVE Data/GISD.csv")
gisd_species <- as.data.frame(gisd$Species)
colnames(gisd_species) <- "gisd"

# our names
names <- unique(eva[, c(2:6)])

# convert species to our name classification
# name
gisd_species$eva <- names$name[match(vegdata::taxname.abbr(gisd_species$gisd),vegdata::taxname.abbr(gsub(" aggr\\.", "", names$name)))]
sum(is.na(gisd_species$eva))
gisd_species$eva[is.na(gisd_species$eva)] <- names$name[match(vegdata::taxname.simplify(gisd_species$gisd[is.na(gisd_species$eva)]),vegdata::taxname.simplify(gsub(" aggr\\.", "", names$name)))]
sum(is.na(gisd_species$eva))
gisd_species$eva[is.na(gisd_species$eva)] <- names$name[match(vegdata::taxname.abbr(gisd_species$gisd[is.na(gisd_species$eva)]),vegdata::taxname.abbr(gsub(" aggr\\.", "", names$species)))]
sum(is.na(gisd_species$eva))

gisd$eva <- gisd_species$eva
#write_csv(gisd,"../EIVE Data/GISD.csv")

impact<- read.csv("I:/Impact_1980.csv")
general <- impact |> 
  group_by(taxa, Neophyte) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots))

general_native <- general[general$Neophyte=="native",]
general <- general[!general$Neophyte=="native",]

gisd$impact <- general$impact[match(gisd$eva, general$taxa)]
gisd$Neophyte <- general$Neophyte[match(gisd$eva, general$taxa)]

gisd$impact[is.na(gisd$impact)]<- general_native$impact[match(gisd$eva[is.na(gisd$impact)], general_native$taxa)]
gisd$Neophyte[is.na(gisd$Neophyte)]<- general_native$Neophyte[match(gisd$eva[is.na(gisd$Neophyte)], general_native$taxa)]

# here we can try whether our species are among the most impactful
Daisy <- read.csv("../GloNAF_Shapefile/Daisy_list.csv", sep=";")
Daisy$taxa <- paste(Daisy$genus, Daisy$specificEpithet)

impact<- read.csv("I:/Impact_1980.csv")
phylo <- read.csv("phylo.csv")
impact[,14:15] <- phylo[match(impact$taxa, phylo$name), 6:7]

check <-vegdata::parse.taxa(unique(impact$taxa))
genus <- check[is.na(check$epi1),]
impact <- impact[!impact$taxa %in% genus$original,]



Daisy$name <- names$name[match(Daisy$taxa, names$name)]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], names$`Matched concept`)]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], names$irena)]
sum(!is.na(Daisy$name))

general <- impact |> 
  group_by(taxa, Neophyte, genus, family) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots))

general_native <- general[general$Neophyte=="native",]
general <- general[!general$Neophyte=="native",]

Daisy$impact <- general$impact[match(Daisy$name, general$taxa)]

Daisy$Neophyte <- general$Neophyte[match(Daisy$name, general$taxa)]

Daisy$impact[is.na(Daisy$impact)] <- general_native$impact[match(Daisy$name[is.na(Daisy$impact)], general_native$taxa)]
Daisy$Neophyte[!is.na(Daisy$impact) & is.na(Daisy$Neophyte)] <- "native"

Daisy$impact[is.na(Daisy$impact)] <- general$impact[match(Daisy$taxa[is.na(Daisy$impact)], general$taxa)]
Daisy$Neophyte[is.na(Daisy$impact)] <- general$Neophyte[match(Daisy$taxa[is.na(Daisy$impact)], general$taxa)]

Daisy$impact[is.na(Daisy$impact)] <- general_native$impact[match(Daisy$taxa[is.na(Daisy$impact)], general_native$taxa)]
Daisy$Neophyte[!is.na(Daisy$impact) & is.na(Daisy$Neophyte)] <- "native"

Daisy |> group_by(Neophyte) |> summarise(n=n(), rel= mean(impact))

Daisy[!duplicated(Daisy$name),] |> group_by(Neophyte) |> summarise(n=n(), rel= mean(impact))

Daisy <- Daisy[!is.na(Daisy$impact),]
Daisy$GISD <- ifelse(Daisy$name %in% gisd$eva, T, F)
colnames(Daisy)
Daisy <- Daisy[, c(20, 21,22,23,17, 9,10)]

Daisy$impact_dom <- impact$RelDiff[match(Daisy$name, impact$taxa[impact$class=="70%-100%"])]

Daisy[!duplicated(Daisy$name),] |> group_by(Neophyte) |> summarise(n=n(), rel= mean(impact), n_high = sum(!is.na(impact_dom)) ,rel_high= mean(impact_dom, na.rm=T))

write_csv(Daisy, "Daisy_impact.csv")
#### 3 Summary ####
# data should be complete from now on, we calculate the number of species in our dataset with eive and div values
colnames(eva)
sum(!is.na(eva$eive_name))/ (length(eva$eive_name))
length(unique(eva$name[!is.na(eva$eive_name)])) / (length(unique(eva$name)))

sum(!is.na(eva$div_name))/ length(eva$div_name)
length(unique(eva$name[!is.na(eva$div_name)])) / (length(unique(eva$name)))
