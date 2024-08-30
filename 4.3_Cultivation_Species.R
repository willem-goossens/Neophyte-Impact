
rm(list=ls())
library(readr)
library(dplyr)

# Load eva and plot data
eva <- read_csv("fullPlotEva_ESy.csv", show_col_types = FALSE)
fullPlotData <- read_csv("fullPlotData_ESy.csv", show_col_types = FALSE)

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

eva$name[eva$name=="Medicago sativa"] <- "Medicago sativa aggr."

# create our own vector
species_vector <- c("Beta vulgaris subsp. vulgaris","Beta vulgaris", "Secale cereale","Triticum aestivum","Pisum sativum", "Solanum tuberosum",
                    "Hordeum vulgare", "Avena sativa", "Zea mays", "Brassica napus","Oryza sativa","Triticum turgidum","Glycine max",
                    "Brassica oleracea", "Sorghum bicolor","Helianthus annuus","Arachis hypogaea","Raphanus sativus", "Medicago sativa aggr.",
                    "Trifolium incarnatum","Phacelia tanacetifolia","Fagopyrum esculentum")

# check number of observations higher than 50
x <- eva[eva$name %in% species_vector & eva$`Cover %`>=50,]
y<- x |> group_by(name) |> summarise(n=n())

# remove plots from both header and eva
plots_cultivated <- eva$PlotObservationID[eva$name %in% species_vector & eva$`Cover %`>=50]
fullPlotData<- fullPlotData[!fullPlotData$PlotObservationID %in% plots_cultivated, ]                    
eva <- eva[!eva$PlotObservationID %in% plots_cultivated,]

# check status species
species_country_status<- read_csv("country_species_ESy.csv", show_col_types = FALSE)
species_vector <- as.data.frame(species_vector)
colnames(species_vector)<- "species"
species_vector$status <- species_country_status$Neophyte[match(species_vector$species, species_country_status$name)]

#write_csv(eva, "fullPlotEva_ESy.csv")
#write_csv(fullPlotData, "fullPlotData_ESy.csv")



#### 2 Invasives ####

gisd <- read.csv("../EIVE Data/GISD.csv", sep = ";")
gisd_species <- gisd$Species

gisd_species[gisd_species %in% uniqueSpecies]

test <- read.csv("coverClassImpactForCandidates_cover_all_layer_cleaned.csv")

test$invasive <- test$taxa %in% gisd$Species

#### 3 Summary ####
# data should be complete from now on, we calculate the number of species in our dataset with eive and div values
eive <- eva[, c(5,8:12)]
sum(!is.na(eva$eive_name))/ (length(eva$eive_name))
length(unique(eva$name[!is.na(eva$eive_name)])) / (length(unique(eva$name))-1)

sum(!is.na(eva$div_name))/ length(eva$div_name)
length(unique(eva$name[!is.na(eva$div_name)])) / (length(unique(eva$name))-1)

length(unique(eva$name))

# we also have to recalculate the total cover for species names seperately, hereby paying attention to the layer in which they are present
