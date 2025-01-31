rm(list=ls())

library(readr)
library(dplyr)
library(ggplot2)
library(stats)

##### 1 READ #####
eva <- read_csv("fullPlotEva_ESy_1980.csv")
fullPlotData <- read_csv("fullPlotData_ESy_1980.csv")

# Data on which species are neophytes
native_intra_analysis=F
if(native_intra_analysis){
  species_country_status<- read_csv("country_species_ESy_1980.csv", show_col_types = FALSE)
} else{
  species_country_status<- read_csv("country_species_ESy_1980.csv", show_col_types = FALSE)
  species_country_status$Neophyte[species_country_status$Neophyte=="native_intra"] <- "native"
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
eva<- left_join(eva, species_country_status[, -c(2:4,6:8)], by= c("Region"= "Region", "name"= "name"))
eva <- eva[!eva$name=="Plant",]


##### 2 DIV #####
###### 2.1 Disturbance Severity #####
test <- aov(Disturbance.Severity ~Neophyte, eva)
summary(test)
TukeyHSD(test)

ggplot(eva, aes(x=Neophyte, y=Disturbance.Severity, color= Neophyte))+
  geom_boxplot()+
  theme_bw()+
  stat_summary(fun= "mean",
              geom = "point", aes(group= Neophyte), size=3)+
  scale_colour_manual(values=c("#1E88E5","#B71C1C", "#FFC107", "#388E3C"), 
                      name="Legend",
                      breaks=c("native","native_intra", "intra","extra"),
                      labels=c("native species","native species alien elsewhere", "intra European neophyte", "extra European neophyte"))


###### 2.2 N #####
test <- aov(EIVEresN ~Neophyte, eva)
summary(test)
TukeyHSD(test)

ggplot(eva, aes(x=Neophyte, y=EIVEresN, color= Neophyte))+
  geom_boxplot()+
  theme_bw()+
  stat_summary(fun= "mean",
               geom = "point", aes(group= Neophyte), size=3)+
  scale_colour_manual(values=c("#1E88E5","#B71C1C", "#FFC107", "#388E3C"), 
                      name="Legend",
                      breaks=c("native","native_intra", "intra","extra"),
                      labels=c("native species","native species alien elsewhere", "intra European neophyte", "extra European neophyte"))



###### 2.3 N niche width #####
test <- aov(EIVEnwN ~Neophyte, eva)
summary(test)
TukeyHSD(test)

ggplot(eva, aes(x=Neophyte, y=EIVEnwN, color= Neophyte))+
  geom_boxplot()+
  theme_bw()+
  stat_summary(fun= "mean",
               geom = "point", aes(group= Neophyte), size=3)+
  scale_colour_manual(values=c("#1E88E5","#B71C1C", "#FFC107", "#388E3C"), 
                      name="Legend",
                      breaks=c("native","native_intra", "intra","extra"),
                      labels=c("native species","native species alien elsewhere", "intra European neophyte", "extra European neophyte"))


##### 3 Databases #####
databases <-  readxl::read_excel("../EVA Data/171_NeophyteInvasions20230216_metadata.xlsx")

data<- fullPlotData |> group_by(Dataset) |> summarise(n=n()) |> mutate(rel_total_in_dataset_used= round(n/length(fullPlotData$PlotObservationID),3))

colnames(databases)[1]<- "Dataset"

databases <- left_join(databases, data)

databases$`# of plots`
databases <- databases |> mutate(rel_used_vs_obtained = round(n/`# of plots`, 3))

#write.csv(databases,"databases.csv", row.names = F)


##### 4. TIME ####
fullPlotData$Date <-  as.Date(fullPlotData$Date, format = "%d.%m.%Y")

# save histogram in jpeg format in current directory
jpeg(file="time histogram.jpeg")

# a histogram we want to save
hist(fullPlotData$Date, breaks= 100, freq=F)

# call this function to save the file 
dev.off()

summary(fullPlotData$Date)

# extra analysis
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


#### 5. Numbers ####
##### 5.1 Species #####
eva_names <- eva |> group_by(name) |> summarise(n =n())
sum(eva_names$n >=30)

eva_names$div <- eva$div_name[match(eva_names$name, eva$name)]
eva_names$eive <- eva$eive_name[match(eva_names$name, eva$name)]

sum(!is.na(eva_names$div))/ nrow(eva_names)
sum(!is.na(eva_names$eive)) / nrow(eva_names)

sum(!is.na(eva$eive_name)) / nrow(eva)
sum(!is.na(eva$div_name)) / nrow(eva)

