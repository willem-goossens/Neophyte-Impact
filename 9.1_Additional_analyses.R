eva <- read_csv("fullPlotEva.csv")
fullPlotData <- read_csv("fullPlotData2.csv")

native_intra_analysis=T
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


# Only observation ID and Region
fullPlot2<- fullPlotData[,c("PlotObservationID","Region")]
# Right join to keep only species present in fullplot (otherwise a lot of NAs)
eva<- right_join(eva, fullPlot2, by = c("PlotObservationID"="PlotObservationID"))
rm(list=c("fullPlot2"))
# Join eva and classification
eva<- left_join(eva, species_country_status, by= c("Region"= "Region", "species"= "species"))

ggplot(eva, aes(x=Disturbance.Severity, fill= Neophyte))+
  geom_histogram(bins=20, position = "dodge", alpha = 0.5,
                 # Note that y = stat(count) is the default behaviour
                 mapping = aes(y = stat(ncount)))

ggplot(eva, aes(x=Neophyte, y=Disturbance.Severity, color= Neophyte))+
  geom_violin()+
  stat_summary(fun= "mean",
              geom = "point", aes(group= Neophyte), size=3)+
  scale_colour_manual(values=c("#1E88E5","#B71C1C", "#FFC107", "#388E3C"), 
                      name="Legend",
                      breaks=c("native","native_intra", "intra","extra"),
                      labels=c("native species","native species alien elsewhere", "intra European neophyte", "extra European neophyte"))



result<-(aov((Disturbance.Severity) ~ Neophyte, eva))
summary(result)
TukeyHSD(result)



databases <-  readxl::read_excel("../EVA Data/171_NeophyteInvasions20230216_metadata.xlsx")

fullPlotData<- read_csv("fullPlotData_cover_all_layer_cleaned.csv")

data<- fullPlotData |> group_by(Dataset) |> summarise(n=n()) |> mutate(rel_total_in_dataset_used= round(n/length(fullPlotData$PlotObservationID),3))

colnames(databases)[1]<- "Dataset"

databases <- left_join(databases, data)

databases$`# of plots`
databases <- databases |> mutate(rel_used_vs_obtained = round(n/`# of plots`, 3))

#write.csv(databases,"databases.csv", row.names = F)

fullPlotData$Date <-  as.Date(fullPlotData$Date, format = "%d.%m.%Y")


# save histogram in jpeg format in current directory
jpeg(file="time histogram.jpeg")

# a histogram we want to save
hist(fullPlotData$Date, breaks= 100, freq=F)

# call this function to save the file 
dev.off()


summary(fullPlotData$Date)

answers_by_date <- fullPlotData %>%
  select(Date, PlotObservationID) %>%
  pad(start_val = min(fullPlotData$Date), end_val = max(fullPlotData$Date), interval = "day") %>%
  group_by(Date) %>%
  summarize(Freq=sum(PlotObservationID)) %>%
  mutate(cumulative = cumsum(replace_na(Freq, 0))) |>mutate(rel_cumsum = cumulative/max(cumulative))

p <- ggplot(answers_by_date, aes(x= Date)) +
  geom_line(aes(y=rel_cumsum), color = "blue")+ scale_x_date(date_breaks = "10 years", date_labels = "%Y")
p

#ggsave("Cumulative time distribution.jpeg", p, width = 25, height = 15, units = "cm")


