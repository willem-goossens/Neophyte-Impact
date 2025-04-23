
library(readr)

# read data
eva <- read_csv("../EVA data/fullPlotEva_new.csv")
fullPlotData <- read_csv("../EVA data/fullPlotData_new.csv")
old <- read.csv("../EVA data/fullPlotData_EUNIS_1980.csv")

fullPlotData <- fullPlotData[!duplicated(fullPlotData$PlotObservationID),]

# downscale model
fast=F
if(fast){
  fullPlotData <- fullPlotData[runif(fullPlotData$PlotObservationID)>0.9,]
  eva <- eva[eva$PlotObservationID %in% fullPlotData$PlotObservationID,]
}

remaining <- fullPlotData[!fullPlotData$PlotObservationID %in% old$PlotObservationID,]
eva <- eva[eva$PlotObservationID %in% remaining$PlotObservationID,]

# now a last check with ESy again
obs <- eva[, c(1,6,9)]
colnames(obs) <- c("RELEVE_NR","TaxonName","Cover_Perc")
any(is.na(obs$TaxonName))

colnames(remaining)
anyNA(remaining$Region)

print(unique(remaining[, c("Country","Region")]), n= 200)

which(duplicated(remaining$PlotObservationID))


remaining[is.na(remaining$Country),]
# now make header
header <- remaining[, c(1, 40, 4,3,5, 44, 43,41, 6)]
header$DUNE[is.na(header$DUNE)]<- "N_DUNES"
colnames(remaining)
str(header)

change <- data.frame(old= c("East.Aegean.Islands","Crete","Republic.of.Ireland", "Baleares","Bosnia.Herzegovina","Corsica",
                            "Czech.Republic","Faroes","Luxemburg","Macedonia","Moldavia","Rf.C",
                            "Rf.CS","Rf.E","Rf.K", "Rf.N",
                            "Rf.NW","Rf.S","Sardinia","Serbia+Kosovo","Sicily","Slovakia",
                            "United.Kingdom"),
                     new= c("Greece","Greece","Ireland","Spain","Bosnia-Herzegovina","France",
                            "Czech Republic","Faroe Islands","Luxembourg", "North Macedonia","Moldova","Russian Federation",
                            "Russian Federation","Russian Federation","Russian Federation","Russian Federation",
                            "Russian Federation","Russian Federation","Spain","Serbia","Italy","Slovak Republic",
                            "United Kingdom"))
header$Region[header$Region %in% change$old] <- change$new[match(header$Region[header$Region %in% change$old], change$old)]
sort(unique(header$Region))

colnames(header)<- c("RELEVE_NR", "Altitude..m",  "DEG_LAT", "DEG_LON", "Country", "Coast_EEA", 'Dunes_Bohn', "Ecoreg", "dataset")
header$`Altitude..m`[header$`Altitude..m` <=0 ] <- 0

header$Ecoreg <- as.numeric(header$Ecoreg)
header$Dunes_Bohn <- as.character(header$Dunes_Bohn)
header$Coast_EEA <- as.character(header$Coast_EEA)


mc <- getOption("mc.cores", 1) 

# run function
source("3.5.1_Bruelheide.R")
# add names to our dataset

types <- data.frame(vegtype.formula.names.short, vegtype.formula.names)

plots <- as.data.frame(result.classification)
colnames(plots) <- "short"
plots$long <- types$vegtype.formula.names[match(plots$short, types$vegtype.formula.names.short)]
plots <- cbind(header$RELEVE_NR, plots)
colnames(plots)[1]<- "ID"

plots$shorter <- sub(" .*","",plots$short)
plots$shorter <- ifelse(nchar(plots$shorter) >= 2, substr(plots$shorter, 1, 2), substr(plots$shorter, 1,1))
plots$shortest <- str_sub(plots$short, 1, 1)


convert <- data.frame(short= c("N ","N1","N2","N3",
                        "Q","Q1","Q2","Q3","Q4","Q5","Qb","Qa",
                        "R","R1","R2","R3","R4","R5","R6","R7",
                        "S","S1","S2","S3","S4","S5","S6","S7","S8","S9","Sa","Sb",
                        "T","T1","T2","T3","T4",
                        "V","V1","V2","V3","V4","V5","V6",
                        "A","A2",
                        "C","C1","C2","C3",
                        "H","H2","H3","H5","H6"),
                      long =c("Coastal habitats","Coastal dunes and sandy shores",
                        "Coastal shingle","Rock cliffs, ledges and shores, including the supralittoral",
                        "Wetlands","Raised and blanket bogs","Valley mires, poor fens and transition mires",
                        "Palsa and polygon mires","Base-rich fens and calcareous spring mires","Helophyte beds","Wetlands","Mires",
                        "Grasslands and lands dominated by forbs, mosses or lichens","Dry grasslands","Mesic grasslands",
                        "Seasonally wet and wet grasslands","Alpine and subalpine grasslands","Woodland fringes and clearings and tall forb stands",
                        "Inland salt steppes","Sparsely wooded grasslands",
                        "Heathlands, scrub and tundra","Tundra","Arctic, alpine and subalpine scrub","Temperate and Mediterranean montane scrub",
                        "Temperate heathland","Maquis, arborescent matorral and thermo-Mediterranean scrub","Garrigue",
                        "Spiny Mediterranean heaths","Thermo-Atlantic xerophytic scrub","Riverine and fen scrub",
                        "Scrub","Dwarf-shrub vegetation",
                        "Forests and other wooded land","Broadleaved deciduous forests","Broadleaved evergreen forests",
                        "Coniferous forests","Lines of trees, small anthropogenic forests, recently felled forest, early-stage forest and coppice",
                        "Vegetated man-made habitats","Arable land and market gardens","*Cultivated areas of gardens and parks",
                        "Artificial grasslands and herb-dominated habitats","*Hedgerows","*Shrub plantations","*Tree dominated man-made habitats",
                        "Marine habitats","Littoral sediment",
                        "Inland surface waters","Surface standing waters","Surface running waters","Littoral zone of inland surface waterbodies",
                        "Inland unvegetated or sparsely vegetated habitats","Screes","Inland cliffs, rock, pavement and outcrops",
                        "Miscellaneous inland habitats with very sparse or no vegetation","Recent volcanic features"
                        ))

plots$type <- convert$long[match(plots$shorter, convert$short)]
plots$habitat <- convert$long[match(plots$shortest, convert$short)]


remaining$habitat <- plots$habitat[match(remaining$PlotObservationID, plots$ID)]

fullPlotData$habitat[!fullPlotData$PlotObservationID %in% old$PlotObservationID] <- remaining$habitat

#write_csv(fullPlotData, "../EVA data/fullPlotData_new.csv")
