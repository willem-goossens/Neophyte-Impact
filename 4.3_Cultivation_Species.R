
rm(list=ls())
library(tidyverse)

# Load eva and plot data
eva <- read_csv("../EVA data/fullPlotEva_4.csv", show_col_types = FALSE)
fullPlotData <- read_csv("../EVA data/fullPlotData_4.csv", show_col_types = FALSE)
cult <- readxl::read_excel("../Extra data/Species names/Cultivated species.xlsx")
species_country_status<- read_csv("../EVA data/country_species_new.csv", show_col_types = FALSE)

# assign species the correct origin
species <- species_country_status[species_country_status$name %in% cult$Species,]
species <- species[!(duplicated(species[,c(1,5,8)])),]
species <- species[, c(1,5,8)]

species$cult <- cult$`Naturalized in non-native range`[match(species$name, cult$Species)]
species$cult[is.na(species$cult)] <- 0

# Downsample by a factor of 100 if fast is selected
fast <- F
if(fast) {
  fullPlotData <- fullPlotData[runif(length(fullPlotData$PlotObservationID)) > 0.90,]
  eva <- eva[eva$PlotObservationID %in% fullPlotData$PlotObservationID,]
}


#### 1 Cultivated ####
# this has to be changed if further

# change crop names
# change crops 
crops <- eva[grepl("crop ", eva$name, ignore.case = TRUE),]

# calculated cover 
cover_total <- eva |> group_by(PlotObservationID) |> summarise(n = sum(`Cover %`))


# get all cultivated species native range
crops <- eva[eva$name %in%  species$name[(species$Neophyte == "native" | species$Neophyte == "native_intra")],]
crops$rel <- crops$`Cover %`/cover_total$n[match(crops$PlotObservationID, cover_total$PlotObservationID)]

x <- crops[crops$rel >= 0.8,]
#print(eva[eva$PlotObservationID == "11590",], n =50)
y<- x |> group_by(name) |> summarise(n=n())


# remove
fullPlotData<- fullPlotData[!fullPlotData$PlotObservationID %in% x$PlotObservationID, ]                    
eva <- eva[!eva$PlotObservationID %in% x$PlotObservationID,]


# get all cultivated species not-native but intra-European
crops <- eva[eva$name %in% species$name[species$Neophyte == "intra"],]
crops$rel <- crops$`Cover %`/cover_total$n[match(crops$PlotObservationID, cover_total$PlotObservationID)]

x <- crops[crops$rel >=0.6,]
y<- x |> group_by(name) |> summarise(n=n())

# remove
fullPlotData<- fullPlotData[!fullPlotData$PlotObservationID %in% x$PlotObservationID, ]  
eva <- eva[!eva$PlotObservationID %in% x$PlotObservationID,]


# get all cultivated species not-native extra-Europea
crops <- eva[eva$name %in% species$name[species$cult=="0" & species$Neophyte == "extra"],]
crops$rel <- crops$`Cover %`/cover_total$n[match(crops$PlotObservationID, cover_total$PlotObservationID)]

x <- crops[crops$rel >=0.4,]
y<- x |> group_by(name) |> summarise(n=n())

# remove
fullPlotData<- fullPlotData[!fullPlotData$PlotObservationID %in% x$PlotObservationID, ]  
eva <- eva[!eva$PlotObservationID %in% x$PlotObservationID,]



# get all cultivated species not-native extra-Europea
crops <- eva[eva$name %in% species$name[species$cult=="1" & species$Neophyte == "extra"],]
crops$rel <- crops$`Cover %`/cover_total$n[match(crops$PlotObservationID, cover_total$PlotObservationID)]

x <- crops[crops$rel >=0.6,]
y<- x |> group_by(name) |> summarise(n=n())

# remove
fullPlotData<- fullPlotData[!fullPlotData$PlotObservationID %in% x$PlotObservationID, ]  
eva <- eva[!eva$PlotObservationID %in% x$PlotObservationID,]



fullPlotData$Date <- as.Date(fullPlotData$Date, format = c("%Y-%m-%d"))
fullPlotData <- fullPlotData[fullPlotData$Date >= "1980-01-01",]

eva <- eva[eva$PlotObservationID %in% fullPlotData$PlotObservationID,]


#### 2 Changes ####
remove <- c("Entodon concinnus","Hygrohypnum luridum")

which(eva$name %in% remove)

eva <- eva[!eva$name %in% remove,]

write_csv(eva, "../EVA data/fullPlotEva_new.csv")
write_csv(fullPlotData, "../EVA data/fullPlotData_new.csv")
#write_csv(y, "../Extra data/Intermediate/removed_cultivated.csv")

#### 3 Time ####
# here again check time
hist(fullPlotData$Date, breaks= 20)
summary(fullPlotData$Date)
sum(fullPlotData$Date < as.Date("1970-01-01"))/ nrow(fullPlotData)

#### 4 Invasives ####
##### 4.1 GISD #####
gisd <- read.csv("../Extra Data/Intermediate/GISD.csv")
gisd_species <- as.data.frame(gisd$Species)
colnames(gisd_species) <- "gisd"

# our names
names <- unique(eva[, c(2:6)])

# convert species to our name classification
# name
gisd_species$eva <- names$name[match((gisd_species$gisd),vegdata::taxname.abbr(gsub(" aggr\\.", "", names$name)))]
sum(is.na(gisd_species$eva))
gisd_species$eva[is.na(gisd_species$eva)] <- names$name[match(vegdata::taxname.simplify(gisd_species$gisd[is.na(gisd_species$eva)]),vegdata::taxname.simplify(gsub(" aggr\\.", "", names$name)))]
sum(is.na(gisd_species$eva))
gisd_species$eva[is.na(gisd_species$eva)] <- names$name[match(vegdata::taxname.abbr(gisd_species$gisd[is.na(gisd_species$eva)]),vegdata::taxname.abbr(gsub(" aggr\\.", "", names$species)))]
sum(is.na(gisd_species$eva))

gisd$eva <- gisd_species$eva
#write_csv(gisd,"../EIVE Data/GISD.csv")

# load impact data
impact<- read.csv("I:/Impact_1980.csv")

check <-vegdata::parse.taxa(unique(impact$taxa))
genus <- check[is.na(check$epi1),]
impact <- impact[!impact$taxa %in% genus$original,]

general <- impact |> 
  group_by(taxa, Neophyte) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots))

general_native <- general[general$Neophyte=="native",]
general <- general[!general$Neophyte=="native",]

gisd$impact <- general$impact[match(gisd$eva, general$taxa)]
gisd$Neophyte <- general$Neophyte[match(gisd$eva, general$taxa)]

gisd$impact[is.na(gisd$impact)]<- general_native$impact[match(gisd$eva[is.na(gisd$impact)], general_native$taxa)]
gisd$Neophyte[is.na(gisd$Neophyte)]<- general_native$Neophyte[match(gisd$eva[is.na(gisd$Neophyte)], general_native$taxa)]


##### 4.2 DAISY ####
# here we can try whether our species are among the most impactful
Daisy <- read.csv("../Extra data/Intermediate/Daisy_list.csv", sep=";")
Daisy$taxa <- paste(Daisy$genus, Daisy$specificEpithet)

# load impact and phylogenetic data
impact<- read.csv("I:/Impact_1980_new.csv")
phylo <- read.csv("../Extra data/Species names/phylo.csv")
impact[,14:15] <- phylo[match(impact$taxa, phylo$name), 6:7]

# remove only genus level
check <-vegdata::parse.taxa(unique(impact$taxa))
genus <- check[is.na(check$epi1),]
impact <- impact[!impact$taxa %in% genus$original,]

# our names
colnames(eva)
names <- eva |> group_by(name, Matched.concept, irena) |> summarise(n=n())

# check names
Daisy$name <- names$name[match(Daisy$taxa, names$name)]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], gsub(" aggr\\.", "", names$name))]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], vegdata::taxname.abbr((names$name)))]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], vegdata::taxname.simplify((names$name)))]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], names$Matched.concept)]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], names$irena)]
sum(!is.na(Daisy$name))

# calculate weighted impact
general <- impact |> 
  group_by(taxa, Neophyte, genus, family) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots))

# get native and alien species seperately
general_native <- general[general$Neophyte=="native",]
general <- general[!general$Neophyte=="native",]

# only assign impact values to alien species in DAISY database
Daisy$impact <- general$impact[match(Daisy$name, general$taxa)]

# assign alien status
Daisy$Neophyte <- general$Neophyte[match(Daisy$name, general$taxa)]

# if there are some native species for which no alien species were observed these can also be added
Daisy$impact[is.na(Daisy$impact)] <- general_native$impact[match(Daisy$name[is.na(Daisy$impact)], general_native$taxa)]
Daisy$Neophyte[!is.na(Daisy$impact) & is.na(Daisy$Neophyte)] <- "native"

# assess relative values
Daisy |> group_by(Neophyte) |> summarise(n=n(), rel= mean(impact))

# sometimes names were duplicated so remove
Daisy[!duplicated(Daisy$name),] |> group_by(Neophyte) |> summarise(n=n(), rel= mean(impact))


# remove non-found species
Daisy <- Daisy[!is.na(Daisy$impact),]
Daisy$GISD <- ifelse(Daisy$name %in% gisd$eva, T, F)
colnames(Daisy)
Daisy <- Daisy[, c(20, 21,22,23,17, 8,9)]

# assign largest value
Daisy$impact_dom <- impact$RelDiff[match(Daisy$name, impact$taxa[impact$class=="70%-100%"])]

# check
Daisy[!duplicated(Daisy$name),] |> group_by(Neophyte) |> summarise(n=n(), rel= mean(impact), n_high = sum(!is.na(impact_dom)) ,rel_high= mean(impact_dom, na.rm=T))

# database
Daisy_plot <- Daisy[!duplicated(Daisy$name),]
Daisy_plot <- Daisy[!Daisy$Neophyte=="native",]

mu <- plyr::ddply(Daisy_plot, "Neophyte", summarise, grp.mean=median(impact))

plot_daisie <- ggplot(Daisy_plot, aes(x= impact, fill= Neophyte))+
  geom_histogram(alpha= 0.75)+
  # Thick white line for halo
  geom_rect(data = mu,
            aes(xmin = grp.mean - 0.05, xmax = grp.mean + 0.05,
            ymin = -Inf, ymax = Inf, fill=Neophyte),
            fill = "white", alpha = 0.2)
  
  # Coloured dashed line overlay
  geom_vline(data = mu, aes(xintercept = grp.mean, color = Neophyte),
             size = 1, linetype = "longdash") +
  ggpubr::theme_pubr()+
  xlab("Overall species impact")+
  ylab("Count")+
  scale_colour_manual(values=c( "#004D40", "#FFC107"), labels=  c("Extra-European",  "Intra-European"))+
  scale_fill_manual(values = c( "#004D40", "#FFC107"), labels=  c("Extra-European",  "Intra-European"))+
  scale_y_continuous(limits=c(0,50))+
  theme(legend.position= c(0.8, 0.85))+
  labs(fill= c("Origin"), colour= c("Origin"))

plot_daisie
#ggsave("../Images/Impact_Daisie.png", plot= plot, width = 5, height = 3)

#write_csv(Daisy, "../Extra data/Results/Daisy_impact.csv")

## 4.3 GloNAF
#glonaf_species<- readxl::read_excel("../Extra data/GLONAF/Taxon_x_List_GLONAF.xlsx")
glonaf_regions<- read_csv("../Extra data/GLONAF/Region_GloNAF_vanKleunenetal2018Ecology.csv", show_col_types = FALSE)

glonaf_species <- read_excel("../Extra data/GLONAF/glonaf_taxon_wcvp.xlsx")
glonaf_region <- read_excel("../Extra data/GLONAF/glonaf_region.xlsx")
glonaf_status <- read_excel("../Extra data/GLONAF/glonaf_flora2.xlsx")

glonaf_region <- glonaf_region[glonaf_region$id %in% glonaf_regions$region_id,]

regions_include <- glonaf_regions$region_id[glonaf_regions$tdwg1_name=="Europe"]
regions_include <- c(regions_include, 889, 349)

glonaf_species <- left_join(glonaf_status, glonaf_species ,by= c("taxon_wcvp_id"="id", "family_wcvp"="family_wcvp"))
glonaf_species <- left_join(glonaf_species, glonaf_regions, by = c("region_id"="region_id"))
europe_species <- glonaf_species[glonaf_species$region_id %in% regions_include,]

glonaf_invasives <- europe_species[europe_species$status=="Invasive",]

eva <- read_csv("../EVA data/fullPlotEVA_new.csv", show_col_types=F)
fullPlotData <- read_csv("../EVA data/fullPlotData_new.csv", show_col_types = FALSE)
species_country_status<- read_csv("../EVA data/country_species_new.csv", show_col_types = FALSE)
species_country_status$Neophyte[species_country_status$Neophyte=="native_intra"] <- "native"

# Only observation ID and Region
fullPlot2<- fullPlotData[,c("PlotObservationID","Region")]
# Right join to keep only species present in fullplot (otherwise a lot of NAs)
eva<- right_join(eva, fullPlot2, by = c("PlotObservationID"="PlotObservationID"))
# Join eva and classification
eva <- left_join(eva, species_country_status[, -c(2:4,6:7)], by= c("Region"= "Region", "name"= "name"))

eva_names <- unique(eva[, c(6,26,25)])
eva_names$glonaf_name <- glonaf_species$taxon_corrected[match(gsub(" aggr.","",eva_names$name), glonaf_species$taxon_corrected)]
colnames(glonaf_invasives)
eva_names <- left_join(eva_names, glonaf_invasives[, c(6,7,12,23)], by = c("glonaf_name"="taxon_corrected","Region"="name"))

#### 5 Summary ####
# data should be complete from now on, we calculate the number of species in our dataset with eive and div values
colnames(eva)
sum(!is.na(eva$eive_name))/ (length(eva$eive_name))
length(unique(eva$name[!is.na(eva$eive_name)])) / (length(unique(eva$name)))

sum(!is.na(eva$div_name))/ length(eva$div_name)
length(unique(eva$name[!is.na(eva$div_name)])) / (length(unique(eva$name)))
