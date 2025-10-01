
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

#write_csv(eva, "../EVA data/fullPlotEva_new.csv")
#write_csv(fullPlotData, "../EVA data/fullPlotData_new.csv")
#write_csv(y, "../Extra data/Intermediate/removed_cultivated.csv")

#### 3 Time ####
# here again check time
hist(fullPlotData$Date, breaks= 20)
summary(fullPlotData$Date)
sum(fullPlotData$Date < as.Date("1970-01-01"))/ nrow(fullPlotData)



#### 4 Invasives ####
library(dplyr)
library(lme4)
library(emmeans)
library(multcomp)
library(ggpubr)
library(ggridges)
library(ggh4x)  
library(tidyverse)

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
impact<- read.csv("I:/Impact_1980_final.csv")

check <-vegdata::parse.taxa(unique(impact$taxa))
genus <- check[is.na(check$epi1),]
impact <- impact[!impact$taxa %in% genus$original,]

general <- impact |> 
  group_by(taxa, Neophyte) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots), n= sum(n))

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
impact<- read.csv("I:/Impact_1980_final.csv")
phylo <- read.csv("../Extra data/Species names/phylo.csv")
impact[,14:15] <- phylo[match(impact$taxa, phylo$name), 6:7]

# remove only genus level
check <-vegdata::parse.taxa(unique(impact$taxa))
genus <- check[is.na(check$epi1),]
impact <- impact[!impact$taxa %in% genus$original,]

Daisy$taxa <- gsub("Fallopia","Reynoutria", Daisy$taxa)

# our names
colnames(eva)
names <- eva |> group_by(name, `Matched concept`, irena) |> summarise(n=n())

library(RSQLite)
## EURO MED
# connect to the database
con <- dbConnect(drv=RSQLite::SQLite(), dbname="C:/Users/u0166342/Documents/Doctoraat/Impact/Extra Data/EURO+MED/EuroSL.sqlite")
# list all tables
tables <- dbListTables(con)
lDataFrames <- vector("list", length=length(tables))
for (i in seq(along=tables)) {
  lDataFrames[[i]] <- dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[i]], "'", sep=""))
}
# retrieve data with euro med classification
Euro <- lDataFrames[[1]]
dbDisconnect(con)
# remove unnecessary data
rm(lDataFrames,tables,con)


# add other data from euro+med
Daisy$taxa2 <- Euro$TaxonConcept[match(vegdata::taxname.abbr(Daisy$taxa), vegdata::taxname.abbr(Euro$TaxonName))]

# create dataframe with species names to feed to the database and to aggregate some species
obs <- Daisy[,c(1,19,2)]
# give necessary names
colnames(obs) <- c("RELEVE_NR","TaxonName","Cover_Perc")
# run function
source("2.2.1_Bruelheide.R")
# add names to our dataset
Daisy$taxa <- obs$TaxonName

# chectaxa# check names
Daisy$name <- names$name[match(Daisy$taxa, names$name)]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], gsub(" aggr\\.", "", names$name))]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], vegdata::taxname.abbr((names$name)))]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], vegdata::taxname.simplify((names$name)))]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], names$`Matched concept`)]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], names$irena)]
sum(!is.na(Daisy$name))
Daisy$name[is.na(Daisy$name)] <- names$name[match(Daisy$taxa[is.na(Daisy$name)], gsub("×", "x", names$name))]
sum(!is.na(Daisy$name))

Daisy <- Daisy[!is.na(Daisy$name),]

# calculate weighted impact
general <- impact |> 
  group_by(taxa, Neophyte, genus, family) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots))

# get native and alien species seperately
general_native <- general[general$Neophyte=="native",]
general <- general[!general$Neophyte=="native",]

# only assign impact values to alien species in DAISY database
Daisy$impact <- general$impact[match(Daisy$name, general$taxa)]
Daisy$impact[is.na(Daisy$impact)] <- general$impact[match(gsub("×", "x", Daisy$name[is.na(Daisy$impact)]), gsub("×", "x", general$taxa))]

# assign alien status
Daisy$Neophyte <- general$Neophyte[match(Daisy$name, general$taxa)]
Daisy$Neophyte[is.na(Daisy$Neophyte)] <- general$Neophyte[match(gsub("×", "x", Daisy$name[is.na(Daisy$Neophyte)]), gsub("×", "x", general$taxa))]


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
Daisy_plot <- Daisy_plot[!Daisy_plot$Neophyte=="native",]

mu <- plyr::ddply(Daisy_plot, "Neophyte", summarise, grp.mean=median(impact))

plot_daisie <- ggplot(Daisy_plot, aes(x= impact, fill= Neophyte))+
  geom_histogram(alpha= 0.75)+
  # Thick white line for halo
  #geom_rect(data = mu,aes(xmin = grp.mean - 0.05, xmax = grp.mean + 0.05, ymin = -Inf, ymax = Inf, fill=Neophyte),fill = "white", alpha = 0.2)
  
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

##### 4.3 GloNAF ####
library(tidyverse)
## ## ##GLONAF DATA  ## ## ##
#glonaf_species<- readxl::read_excel("../Extra data/GLONAF/Taxon_x_List_GLONAF.xlsx")
glonaf_regions<- read_csv("../Extra data/GLONAF/Region_GloNAF_vanKleunenetal2018Ecology.csv", show_col_types = FALSE)

library(readxl)
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

length(unique(glonaf_invasives$name))

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

eva_names <- unique(eva[, c(6,25,24)])
eva_names$name[eva_names$name=="Hippophaë rhamnoides"] <- "Hippophae rhamnoides"
eva_names$glonaf_name <- glonaf_species$taxa_accepted[match(gsub(" aggr.","",eva_names$name), glonaf_species$taxa_accepted)]
eva_names$glonaf_name[is.na(eva_names$glonaf_name)] <- glonaf_species$taxa_accepted[match(gsub("×","x",eva_names$name[is.na(eva_names$glonaf_name)]), glonaf_species$taxa_accepted)]
eva_names$glonaf_name[is.na(eva_names$glonaf_name)] <- glonaf_species$taxa_accepted[match(gsub(" aggr.","",eva_names$name[is.na(eva_names$glonaf_name)]), glonaf_species$taxon_corrected)]
colnames(glonaf_invasives)
glonaf_sum <- glonaf_invasives |> group_by(taxa_accepted, status) |> summarise(n = n() )
glonaf_sum <- glonaf_sum[!is.na(glonaf_sum$taxa_accepted),]

eva_names$invasive <- NA
eva_names$invasive[eva_names$glonaf_name %in% glonaf_sum$taxa_accepted] <- "GloNAF"
eva_names$n <- glonaf_sum$n[match(eva_names$glonaf_name, glonaf_sum$taxa_accepted)]

eva_names_n <- eva_names

# left_join(eva_names, glonaf_sum, by = c("glonaf_name"="taxon_corrected"))

## ## ## IMPACT DATA  ## ## ##
## IMPACT
# load impact and phylogenetic data
impact<- read.csv("I:/Impact_1980_final.csv")
phylo <- read.csv("../Extra data/Species names/phylo.csv")
impact[,14:15] <- phylo[match(impact$taxa, phylo$name), 6:7]
impact$taxa[impact$taxa=="Hippophaë rhamnoides"] <- "Hippophae rhamnoides"

# remove only genus level
check <-vegdata::parse.taxa(unique(impact$taxa))
genus <- check[is.na(check$epi1),]
impact <- impact[!impact$taxa %in% genus$original,]

# our names
colnames(eva)
names <- eva |> group_by(name, `Matched concept`, irena) |> summarise(n=n())



# calculate weighted impact
general <- impact |> 
  group_by(taxa, Neophyte, genus, family) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots))


## ## ## Combine ## ## ##
eva_names <- eva_names[!is.na(eva_names$invasive),]
eva_names <- eva_names[!eva_names$Neophyte=="native",]

general$invasive <- NA
general$invasive[general$Neophyte != "native"] <- eva_names$invasive[match(general$taxa[general$Neophyte != "native"], eva_names$name)]

table(general$Neophyte[(!is.na(general$invasive))])


## ## ## Plot ## ## ##
plot_glonaf <- general[!is.na(general$invasive),]

mu <- plyr::ddply(plot_glonaf, "Neophyte", summarise, grp.mean=median(impact))

plot_daisie <- ggplot(plot_glonaf, aes(x= impact, fill= Neophyte))+
  geom_histogram(alpha= 0.75)+
  # Thick white line for halo
  #geom_rect(data = mu, aes(xmin = grp.mean - 0.05, xmax = grp.mean + 0.05, ymin = -Inf, ymax = Inf, fill = Neophyte), alpha = 0.2)+
  # Coloured dashed line overlay
  geom_vline(data = mu, aes(xintercept = grp.mean, color = Neophyte),
           size = 1, linetype = "longdash") +
  ggpubr::theme_pubr()+
  xlab("Overall species impact")+
  ylab("Count")+
  scale_colour_manual(values=c( "#004D40", "#FFC107"), labels=  c("Extra-European",  "Intra-European"))+
  scale_fill_manual(values = c( "#004D40", "#FFC107"), labels=  c("Extra-European",  "Intra-European"))+
  scale_y_continuous(limits=c(0,13))+
  theme(legend.position= c(0.8, 0.85))+
  labs(fill= c("Origin"), colour= c("Origin"))

plot_daisie

length(Daisy_plot$taxa2)
Daisy_plot

## ## ## IMPACT DATA DAISIE ## ## ##
plot_data <- general %>%
  mutate(
    daisie   = ifelse(taxa %in% Daisy_plot$taxa2 & Neophyte != "native", "Daisie", NA),
    invasive = ifelse(Neophyte == "native", NA, invasive),
    type     = "All"
  ) %>%
  bind_rows(
    filter(general, taxa %in% Daisy_plot$name & Neophyte != "native") %>% mutate(type = "Daisie"),
    filter(general, !is.na(invasive)) %>% mutate(type = "GloNAF")
  ) %>%
  mutate(combined = paste(type, Neophyte))

## quick check
colnames(plot_data)

## ANOVA & posthoc
kw <- aov(impact ~ combined, data = plot_data)
TukeyHSD(kw)

## GLMM
model <- lmer(impact ~ -1 + combined + (1|family/genus), data = plot_data)
summary(model)

emm <- emmeans(model, pairwise ~ combined)
cld_results <- cld(emm$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey", reversed = TRUE) %>%
  mutate(.group = gsub(" ", "", .group))

## prep for plotting
y_max <- max(plot_data$impact, na.rm = TRUE)
y_n   <- y_max * 1.1
y_let <- y_max * 1.2

## sample sizes
n_counts <- plot_data %>%
  group_by(combined) %>%
  summarise(n = n(), .groups = "drop")

## merge letters with levels
letters_df <- cld_results[,c(1,7)]

## Plot
ggplot(plot_data, aes(x = combined, y = impact, fill = Neophyte)) +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.25, alpha = 0.8, fill = "white", outlier.shape = NA) +
  stat_summary(fun = "mean", geom = "point", aes(group = Neophyte), size = 3) +
  theme_pubr() +
  scale_colour_manual(values = c("#004D40", "#FFC107", "#1E88E5", "darkgreen")) +
  scale_fill_manual(values = c("#004D40", "#FFC107", "#1E88E5", "darkgreen")) +
  guides(color = "none", fill = "none") +
  scale_x_discrete(labels = c(
    "All Extra-\nEuropean", "All Intra-\nEuropean", "All Native",
    "DAISIE Extra-\nEuropean", "DAISIE Intra-\nEuropean",
    "GloNAF Extra-\nEuropean", "GloNAF Intra-\nEuropean"
  )) +
  theme(
    axis.text.x  = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 10),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(min(plot_data$impact), y_let * 1.1)) +
  geom_text(
    data = n_counts,
    aes(x = combined, y = y_n, label = paste0("n=", n)),
    inherit.aes = FALSE, size = 4, vjust = 0
  ) +
  geom_text(
    data = letters_df,
    aes(x = combined, y = y_let, label = .group),
    inherit.aes = FALSE, size = 4, vjust = -0.5
  )

 # for nested facets

## --- Data prep ---
plot_data <- general %>%
  mutate(
    daisie   = ifelse(taxa %in% Daisy_plot$name & Neophyte != "native", "Daisie", NA),
    invasive = ifelse(Neophyte == "native", NA, invasive),
    type     = "All"
  ) %>%
  bind_rows(
    filter(general, taxa %in% Daisy_plot$name & Neophyte != "native") %>% mutate(type = "Daisie"),
    filter(general, !is.na(invasive)) %>% mutate(type = "GloNAF")
  ) %>%
  mutate(combined = paste(type, Neophyte))

## --- Model & significance ---
model <- lmer(impact ~ -1 + combined + (1|family/genus), data = plot_data)
emm   <- emmeans(model, pairwise ~ combined)
cld_results <- cld(emm$emmeans, alpha = 0.05, Letters = letters,
                   adjust = "tukey", reversed = TRUE) %>%
  mutate(.group = gsub(" ", "", .group))

letters_df <- cld_results[,c(1,7)] %>%
  rename(combined = combined, letter = .group)

n_counts <- plot_data %>%
  group_by(combined) %>%
  summarise(n = n(), .groups = "drop")

plot_labels <- left_join(letters_df, n_counts, by = "combined") 

## --- Add dummy spacer rows ---
spacers <- tibble(
  impact   = NA,              # no value
  Neophyte = NA,              # no colour
  type     = c("Spacer1", "Spacer2"),
  combined = c("Spacer1", "Spacer2")
)

plot_data <- bind_rows(plot_data, spacers)

plot_data2 <- plot_data %>%
  mutate(combined = factor(combined, levels = rev(c(
    "All extra", "All intra", "All native",
    "Spacer1",
    "Daisie extra", "Daisie intra",
    "Spacer2",
    "GloNAF extra", "GloNAF intra"
  ))))


plot_invasive <- ggplot(plot_data2, aes(x = impact, y = combined, fill = Neophyte)) +
  stat_density_ridges(alpha = 0.5, rel_min_height = 0.01, scale = 1.2,
                      quantile_lines = TRUE, quantile_fun = mean, na.rm = TRUE) +
  theme_ridges() +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  ) +
  scale_fill_manual(values = c("#004D40", "#FFC107", "#1E88E5")) +
  scale_y_discrete(labels = rev(c(
    "All Extra-European", "All Intra-European", "All Native",
    "",   # spacer
    "DAISIE Extra-European", "DAISIE Intra-European",
    "",   # spacer
    "GloNAF Extra-European", "GloNAF Intra-European"
  ))) +
  labs(x = "Impact")+
  geom_text(
    data = plot_labels,
    aes(x = max(plot_data2$impact, na.rm = TRUE) * 0.75, 
        y = combined, label = letter),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  ## sample sizes
  geom_text(
    data = plot_labels,
    aes(x = max(plot_data2$impact, na.rm = TRUE) * 0.70, 
        y = combined, label = paste0("n=", n)),
    inherit.aes = FALSE,
    hjust = 1,
    size = 4
  ) +
  coord_cartesian(xlim = c(min(plot_data$impact, na.rm = TRUE),
                           max(plot_data$impact, na.rm = TRUE) * 0.80))+
  geom_vline(xintercept=0, linetype='dotted', alpha=1)+
  theme(axis.ticks.y = element_blank())

plot_invasive

ggsave(plot=plot_invasive , "../Images/Invasive_status.png", width=8, height=5)


sum(Daisy_plot$impact[Daisy_plot$Neophyte=="intra"]>=0)/length(Daisy_plot$impact[Daisy_plot$Neophyte=="intra"]>=0)
sum(Daisy_plot$impact[Daisy_plot$Neophyte=="intra"]<0)
sum(Daisy_plot$impact[Daisy_plot$Neophyte=="extra"]>=0)/length(Daisy_plot$impact[Daisy_plot$Neophyte=="extra"]>=0)

sum(plot_glonaf$impact[plot_glonaf$Neophyte=="intra"]>=0)/length(plot_glonaf$impact[plot_glonaf$Neophyte=="intra"]>=0)
sum(plot_glonaf$impact[plot_glonaf$Neophyte=="extra"]>=0)/length(plot_glonaf$impact[plot_glonaf$Neophyte=="extra"]>=0)


# create database for Jurgen and co

general$Daisie <- NA
general$Daisie[general$taxa %in% Daisy_plot$name] <- "DAISIE"
colnames(general) <- c("Name","Status","Genus","Family","Overall impact", "GloNAF","DAISIE")
general$GloNAF[!is.na(general$GloNAF)] <- "GloNAF"
general$DAISIE[general$Status=="native"] <- NA

general <- general[, c(1,3,4,2,5,6,7)]

general$Status[general$Status=="native"] <- "Native"
general$Status[general$Status=="intra"] <- "Intra-European alien"
general$Status[general$Status=="extra"] <- "Extra-European alien"
general$Status[general$Status=="Native" & general$Name %in% general$Name[general$Status =="Intra-European alien"]] <- "Native alien elsewhere"

eva_names_n <- eva_names_n[!duplicated(eva_names_n[, c(1)]),]
general$n <- eva_names_n$n[match(general$Name, eva_names_n$name)]
general$n[general$Status=="Native"] <- NA
general$n[general$Status=="Native alien elsewhere"] <- NA


general$GloNAF[!is.na(general$n)] <- paste(general$GloNAF[!is.na(general$n)], paste("(",general$n[!is.na(general$n)], ")", sep = ""))

general <- general[, c(1:7)]

speciesDominance<- read.csv("../Results/speciesDominance_1980.csv")
speciesDominance$neophyte[speciesDominance$names %in% speciesDominance$names[speciesDominance$neophyte=="intra"] & speciesDominance$neophyte=="native"] <- "native_intra"


names <- speciesDominance$names
check <-vegdata::parse.taxa(names)
genus <- check[is.na(check$epi1),]
speciesDominance <- speciesDominance[!speciesDominance$names %in% genus$original,]
native_dom <- speciesDominance[speciesDominance$neophyte=="native",]
intra_dom <- speciesDominance[speciesDominance$neophyte=="intra",]
extra_dom <- speciesDominance[speciesDominance$neophyte=="extra",]
native_intra_dom <- speciesDominance[speciesDominance$neophyte=="native_intra",]

general$Frequency <- NA
general$Dominance <- NA

general[general$Status=="Native", 8:9] <- native_dom[match(general$Name[general$Status=="Native"], native_dom$names), c(5,3)]
general[general$Status=="Intra-European alien", 8:9] <- intra_dom[match(general$Name[general$Status=="Intra-European alien"], intra_dom$names), c(5,3)]
general[general$Status=="Extra-European alien", 8:9] <- extra_dom[match(general$Name[general$Status=="Extra-European alien"], extra_dom$names), c(5,3)]
general[general$Status=="Native alien elsewhere", 8:9] <- native_intra_dom[match(general$Name[general$Status=="Native alien elsewhere"], native_intra_dom$names), c(5,3)]

general[, c(5,9)] <- round(general[, c(5,9)], 2)


#write_csv(general, "../Results/Overall_impact.csv")


#### 5 Summary ####
# data should be complete from now on, we calculate the number of species in our dataset with eive and div values
colnames(eva)
sum(!is.na(eva$eive_name))/ (length(eva$eive_name))
length(unique(eva$name[!is.na(eva$eive_name)])) / (length(unique(eva$name)))

sum(!is.na(eva$div_name))/ length(eva$div_name)
length(unique(eva$name[!is.na(eva$div_name)])) / (length(unique(eva$name)))



#### 6 Native intra ####
plot_data2$Neophyte[plot_data2$Neophyte=="native" & plot_data2$taxa %in% plot_data2$taxa[plot_data2$Neophyte=="intra"]] <- "native_intra"
plot_data2$combined <- as.character(plot_data2$combined)
plot_data2$combined[plot_data2$combined=="All native" & plot_data2$taxa %in% plot_data2$taxa[plot_data2$Neophyte=="intra"]] <- (c("All native intra"))

plot_data2 <- plot_data2 %>%
  mutate(combined = factor(combined, levels = rev(c(
    "All extra", "All intra", "All native", "All native intra",
    "Spacer1",
    "Daisie extra", "Daisie intra",
    "Spacer2",
    "GloNAF extra", "GloNAF intra"
  ))))


model <- lmer(impact ~ -1 + combined + (1|family/genus/taxa), data = plot_data2)
DHARma::
emm   <- emmeans(model, pairwise ~ combined)
cld_results <- cld(emm$emmeans, alpha = 0.05, Letters = letters,
                   adjust = "tukey", reversed = TRUE) %>%
  mutate(.group = gsub(" ", "", .group))

letters_df <- cld_results[,c(1,7)] %>%
  rename(combined = combined, letter = .group)

n_counts <- plot_data2 %>%
  group_by(combined) %>%
  summarise(n = n(), .groups = "drop")

plot_labels <- left_join(letters_df, n_counts, by = "combined") 

plot_data3 <- plot_data2

plot_data2 <- plot_data2[plot_data2$impact < max(plot_data2$impact, na.rm = TRUE) * 0.70 | is.na(plot_data2$impact), ]

plot_invasive <- ggplot(plot_data2, aes(x = impact, y = combined, fill = Neophyte)) +
  stat_density_ridges(alpha = 0.5, rel_min_height = 0.01, scale = 1.2,
                      quantile_lines = TRUE, quantile_fun = mean, na.rm = TRUE) +
  theme_ridges() +
  theme_pubr() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.text.y  = element_text(size = 12),
    axis.text.x  = element_text(size = 12)
  ) +
  scale_fill_manual(values = c("#004D40", "#FFC107", "#1E88E5", "darkgreen")) +
  scale_y_discrete(labels = rev(c(
    "All Extra-European", "All Intra-European", "All native", "All native alien elsewhere",
    "",   # spacer
    "DAISIE Extra-European", "DAISIE Intra-European",
    "",   # spacer
    "GloNAF Extra-European", "GloNAF Intra-European"
  ))) +
  labs(x = "Impact")+
  geom_text(
    data = plot_labels,
    aes(x = max(plot_data3$impact, na.rm = TRUE) * 0.75, 
        y = combined, label = letter),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4
  ) +
  ## sample sizes
  geom_text(
    data = plot_labels,
    aes(x = max(plot_data3$impact, na.rm = TRUE) * 0.70, 
        y = combined, label = paste0("n=", n)),
    inherit.aes = FALSE,
    hjust = 1,
    size = 4
  ) +
  coord_cartesian(xlim = c(min(plot_data$impact, na.rm = TRUE),
                           max(plot_data$impact, na.rm = TRUE) * 0.80))+
  geom_vline(xintercept=0, linetype='dotted', alpha=1)+
  theme(axis.ticks.y = element_blank())

plot_invasive

ggsave(plot=plot_invasive , "../Images/Invasive_status_native_intra.png", width=8, height=5)

test <- general[(is.na(general$GloNAF) & is.na(general$DAISIE)),]

test <- test[!(test$Status=="Native" | test$Status =="Native alien elsewhere"),]

round(mean(test$`Overall impact`), 4)
