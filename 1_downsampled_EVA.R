#Load packages
library(readr)
library(tidyr)

#Load EVA data
eva <- read_delim("../EVA Data/171_NeophyteInvasions20230216_notJUICE_species.csv","\t")

# Function needed later
# Function which allows computing the number of unique plot IDs found in dataframes structured like EVA
nUniquePlotIDs <- function(evaDB) {
  length(unique(evaDB$PlotObservationID))
}

#Random boleans (False or True) vector with propabilities of larger than 0.99 or 0.999
# Generate random booleans indicating which Plot IDs to include
largestPlotID <- eva$PlotObservationID[length(eva$Taxonomy)]
onePercentOfPlotIDs <- runif(largestPlotID) > 0.99
onePromilleOfPlotIDs <- runif(largestPlotID) > 0.999


# Select 1% or 0.1% of the plots from EVA
onePercentEva <- eva[onePercentOfPlotIDs[eva$PlotObservationID],]
onePromilleEva <- eva[onePromilleOfPlotIDs[eva$PlotObservationID],]


# Verify that the expected number of plots was actually selecte
numberOfUniquePlotIDsInEva <- nUniquePlotIDs(eva)
fractionOfPlotsInOnePercentEva <- numberOfUniquePlotIDsInEva / nUniquePlotIDs(onePercentEva)
# Should be about 100
fractionOfPlotsInOnePercentEva
fractionOfPlotsInOnePromilleEva <- numberOfUniquePlotIDsInEva / nUniquePlotIDs(onePromilleEva)
# Should be about 1000
fractionOfPlotsInOnePromilleEva

# Store the resulting data in a file
#write_delim(onePercentEva, "../EVA Data/onePercentEva.csv", delim = "\t")
#write_delim(onePromilleEva, "../EVA Data/onePromilleEva.csv", delim = "\t")

