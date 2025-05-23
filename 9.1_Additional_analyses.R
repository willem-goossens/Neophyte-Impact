rm(list=ls())

library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
##### 1 READ #####
eva <- read_csv("../EVA data/fullPlotEva_new.csv")
fullPlotData <- read_csv("../EVA data/fullPlotData_new.csv")

# Data on which species are neophytes
native_intra_analysis=F
if(native_intra_analysis){
  species_country_status<- read_csv("../EVA data/country_species_new.csv", show_col_types = FALSE)
} else{
  species_country_status<- read_csv("../EVA data/country_species_new.csv", show_col_types = FALSE)
  species_country_status$Neophyte[species_country_status$Neophyte=="native_intra"] <- "native"
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
  padr::pad(start_val = min(fullPlotData$Date), end_val = max(fullPlotData$Date), interval = "day") %>%
  group_by(Date) %>%
  summarize(Freq=sum(PlotObservationID)) %>%
  mutate(cumulative = cumsum(replace_na(Freq, 0))) |>mutate(rel_cumsum = cumulative/max(cumulative))

p <- ggplot(answers_by_date, aes(x= Date)) +
  geom_line(aes(y=rel_cumsum), color = "blue")+ scale_x_date(date_breaks = "10 years", date_labels = "%Y")
p


#### 5. Numbers ####
###### 5.1 Species ######
eva_names <- eva |> group_by(name, Neophyte) |> summarise(n =n())
sum(eva_names$n >=30)

analysed <- eva_names[eva_names$n>=30,]

eva_names$div <- eva$div_name[match(eva_names$name, eva$name)]
eva_names$eive <- eva$eive_name[match(eva_names$name, eva$name)]

sum(!is.na(eva_names$div))/ nrow(eva_names)
sum(!is.na(eva_names$eive))/ nrow(eva_names)

sum(!is.na(eva$eive_name)) / nrow(eva)
sum(!is.na(eva$div_name)) / nrow(eva)

sum(eva_names$Neophyte=="extra")
# read species dominance, which is more correct
species_dominance <- read.csv("../Results/speciesDominance_1980.csv")

setdiff(species_dominance$names, analysed$name)
setdiff( analysed$name, species_dominance$names)

#### 6. Traits #####
###### 6.1 Plot EIVE #######
colnames(eva)

eva2 <- eva |> group_by(name, Neophyte,eive_name, div_name, EIVEresM, EIVEresN, EIVEresL, EIVEresR, EIVEresT, Disturbance.Severity, Disturbance.Frequency,
                        Grazing.Pressure, Mowing.Frequency, Soil.Disturbance) |> summarise(n=n())

eva2[duplicated(eva2[,1:2]) | duplicated(eva2[,1:2], fromLast=T) ,]
table(eva2$Neophyte)

sum(!is.na(eva2$eive_name))/length(eva2$name)
sum(!is.na(eva$eive_name))/length(eva$PlotObservationID)
sum(!is.na(eva2$div_name))/length(eva2$name)
sum(!is.na(eva$div_name))/length(eva$PlotObservationID)

names <- c("EIVE-M","EIVE-N","EIVE-L", "EIVE-R","EIVE-T","Disturbance Severity","Disturbance Frequency","Grazing Pressure", "Mowing Frequency","Soil Disturbance")
phylo <- read_csv("../Extra data/Species names/phylo.csv", show_col_types = F)

library(emmeans)
library(multcomp)

for(i in 5: 14){

var <- colnames(eva2)[i]
eva3 <- eva2[, c(1,2,i)]
colnames(eva3)[3] <- "var"
eva3 <- eva3[!is.na(eva3$var),]

eva3[, c(4:5)] <- phylo[match(eva3$name, phylo$name), 6:7]


# glmm
model <- lmer(var~ -1+ Neophyte + (1|family/genus), data= eva3 )
summary(model)
emm<- emmeans::emmeans(model, pairwise ~ Neophyte)

if(i < 1){
  y_n <- 10
  y_let <- 10.2
} else {
  y_n <- max(eva3$var)+max(eva3$var)/10
  y_let <- max(eva3$var)+max(eva3$var)/5
}

class(emm$emmeans)

cld_results <- cld(emm$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey", reversed=T)
cld_results$.group <- gsub(" ","",cld_results$.group)
cld_results

plot <- ggplot(eva3, aes(x=Neophyte, y= var, fill=Neophyte))+
  geom_violin(alpha=0.5, scale="width") +
  geom_boxplot(width=0.25, alpha=0.8, fill="white", outlier.shape = NA) +
  theme_pubr() +
  stat_summary(fun="mean", geom="point", aes(group=Neophyte), size=3) +
  scale_colour_manual(values=c( "#004D40","#FFC107", "#1E88E5", "darkgreen")) +
  scale_fill_manual(values=c( "#004D40","#FFC107", "#1E88E5", "darkgreen")) +
  guides(color="none") +
  theme(legend.position="none") +
  scale_x_discrete(labels=c("Extra-\nEuro pean", "Intra-\nEuropean", "Native", "Native\nalien elsewhere")) +
  theme(axis.text.x=element_text(size=12), axis.title.y=element_text(size=12), axis.title.x=element_blank(), 
        axis.text.y=element_text(size=10)) +
  scale_y_continuous(limits=c(0, y_let+y_let/10))+
  ylab(names[i-4])+
  annotate("text", x=1, y=y_n, label=paste("n=", length(eva3$name[eva3$Neophyte=="extra"])), size=4, hjust=0.5)+
  annotate("text", x=2, y=y_n, label=paste("n=", length(eva3$name[eva3$Neophyte=="intra"])), size=4, hjust=0.5) +
  annotate("text", x=3, y=y_n, label=paste("n=", length(eva3$name[eva3$Neophyte=="native"])), size=4, hjust=0.5) +
  annotate("text", x = 1, y = y_let, label = cld_results$.group[cld_results$Neophyte=="extra"], size = 4,vjust = -0.5, hjust = 0.5, alpha=1) +
  annotate("text", x = 2, y = y_let, label = cld_results$.group[cld_results$Neophyte=="intra"], size = 4,vjust = -0.5, hjust = 0.5, alpha=1) +
  annotate("text", x = 3, y = y_let, label = cld_results$.group[cld_results$Neophyte=="native"], size = 4,vjust = -0.5, hjust = 0.5, alpha=1) 
 plot

 assign(gsub(" ","_",var), plot) 
 
}

plot <- ggarrange(EIVEresN, EIVEresR, EIVEresL, EIVEresL, EIVEresM,Disturbance.Severity,  Disturbance.Frequency, 
                  Grazing.Pressure, Mowing.Frequency, Soil.Disturbance, nrow=5, ncol=2, labels = LETTERS[1:10],
                  font.label = list(size = 12), align="hv")
plot

#ggsave("../Images/Species_values.svg", plot= plot, width = 10, height = 10)


###### 6.2 Plot Traits #######
eva_names <- eva |> group_by(name, Neophyte) |>summarise(n=n())

Diaz <- read_csv("../Extra data/Traits/Eva_diaz.csv")
eva_names <- left_join(eva_names, Diaz[!duplicated(Diaz$name),-c(2:6)], by=c("name"="name"))

mean <- read_csv("../Results/trait_normalized.csv")

eva_names$LA[is.na(eva_names$LA)] <- mean$LA[match(eva_names$name[is.na(eva_names$LA)], mean$SpeciesName)]
eva_names$LMA[is.na(eva_names$LMA)] <- mean$LMA[match(eva_names$name[is.na(eva_names$LMA)], mean$SpeciesName)]
eva_names$SSD[is.na(eva_names$SSD)] <- mean$SSD[match(eva_names$name[is.na(eva_names$SSD)], mean$SpeciesName)]
eva_names$SM[is.na(eva_names$SM)] <- mean$SM[match(eva_names$name[is.na(eva_names$SM)], mean$SpeciesName)]
eva_names$H[is.na(eva_names$H)] <- mean$H[match(eva_names$name[is.na(eva_names$H)], mean$SpeciesName)]
eva_names$Nmass[is.na(eva_names$Nmass)] <- mean$Nmass[match(eva_names$name[is.na(eva_names$Nmass)], mean$SpeciesName)]

# standardize
test <- eva_names
test[, c(4:9)] <- log(test[,c(4:9)])
test[, c(4:9)] <- scale(test[, c(4:9)])
j=4
# remove outliers (see Diaz et al 2016 for number 4)
for(j in 4:9){
  tr <- as.matrix(test[,j])
  print(sum(is.na(tr)))
  tr <- ifelse(tr <= -4, NA, tr)
  tr <- ifelse(tr >= 4, NA, tr)
  print(sum(is.na(tr)))
  which(is.na(tr))
  eva_names[is.na(tr),j] <- NA
}

phylo <- read.csv("../Extra data/Species names/phylo.csv")
eva_names[, c(14,13)] <- phylo[match(eva_names$name, phylo$name), 6:7]
names <- c("Leaf area [mm²]", "LMA [g/m²]","SSD [mg/mm³]","Seed mass [mg]","Leaf N [mg/g]","Height [m]")



i=4
for(i in 4: 9){
  
  var <- colnames(eva_names)[i]
  eva3 <- eva_names[, c(1,2,i, 14,13)]
  colnames(eva3)[3:5] <- c("var","genus","family")
  eva3 <- eva3[!is.na(eva3$var),]
  
  
  # glmm
  model <- lmer(var~ -1+Neophyte + (1|family/genus), data= eva3 )
  summary(model)
  emm<- emmeans::emmeans(model, pairwise ~ Neophyte)
  
  if(i < 1){
    y_n <- 10
    y_let <- 10.2
  } else {
    value <- as.numeric(quantile(eva3$var)[4] + IQR(eva3$var)*1.5)
    y_n <- value + value/10
    y_let <- value+ value/5
  }

  cld_results <- cld(emm$emmeans, alpha = 0.05, Letters = letters, adjust = "tukey", reversed=T)
  cld_results$.group <- gsub(" ","",cld_results$.group)
  cld_results
  
  plot <- ggplot(eva3, aes(x=Neophyte, y= var, fill=Neophyte))+
    geom_violin(alpha=0.5, scale="width") +
    geom_boxplot(width=0.25, alpha=0.8, fill="white", outlier.shape = NA) +
    theme_pubr() +
    stat_summary(fun="mean", geom="point", aes(group=Neophyte), size=3) +
    scale_colour_manual(values=c( "#004D40","#FFC107", "#1E88E5", "darkgreen")) +
    scale_fill_manual(values=c( "#004D40","#FFC107", "#1E88E5", "darkgreen")) +
    guides(color="none") +
    theme(legend.position="none") +
    scale_x_discrete(labels=c("Extra-\nEuropean", "Intra-\nEuropean", "Native", "Native\nalien elsewhere")) +
    theme(axis.text.x=element_text(size=12), axis.title.y=element_text(size=12), axis.title.x=element_blank(), 
          axis.text.y=element_text(size=10)) +
    scale_y_continuous(limits=c(0, y_let+y_let/10))+
    ylab(names[i-3])+
    annotate("text", x=1, y=y_n, label=paste("n=", length(eva3$taxa[eva3$Neophyte=="extra"])), size=4, hjust=0.5)+
    annotate("text", x=2, y=y_n, label=paste("n=", length(eva3$taxa[eva3$Neophyte=="intra"])), size=4, hjust=0.5) +
    annotate("text", x=3, y=y_n, label=paste("n=", length(eva3$taxa[eva3$Neophyte=="native"])), size=4, hjust=0.5) +
    annotate("text", x = 1, y = y_let, label = cld_results$.group[cld_results$Neophyte=="extra"], size = 4,vjust = -0.5, hjust = 0.5, alpha=1) +
    annotate("text", x = 2, y = y_let, label = cld_results$.group[cld_results$Neophyte=="intra"], size = 4,vjust = -0.5, hjust = 0.5, alpha=1) +
    annotate("text", x = 3, y = y_let, label = cld_results$.group[cld_results$Neophyte=="native"], size = 4,vjust = -0.5, hjust = 0.5, alpha=1) 
  plot
  
  assign(gsub(" ","_",var), plot) 
  
}

plot <- ggarrange(LA,SSD,LMA,SM,Nmass,H, nrow=3, ncol=2, labels = LETTERS[1:10],
                  font.label = list(size = 12), align="hv")
plot

#ggsave("../Images/Species_traits.svg", plot= plot, width = 10, height = 6)


##### 7 Impact #####

x <- read_csv('I:/Impact_1980_new.csv')
speciesDominance<- read.csv("../Results/speciesDominance_1980.csv")

# remove genus data
check <-vegdata::parse.taxa(unique(x$taxa))
genus <- check[is.na(check$epi1),]
x <- x[!x$taxa %in% genus$original,]

phylo <- read_csv("../Extra data/Species names/phylo.csv", show_col_types = F)
speciesDominance[, c(7:8)] <- phylo[match(speciesDominance$names, phylo$name), 6:7]

names <- speciesDominance$names
check <-vegdata::parse.taxa(names)
genus <- check[is.na(check$epi1),]
speciesDominance <- speciesDominance[!speciesDominance$names %in% genus$original,]

general <- x |> 
  group_by(taxa, Neophyte) |> 
  summarise(impact= sum((RelDiff*n)/numberOfPlots))


general$number <- speciesDominance$numberOfPlot[match(general$taxa, speciesDominance$names)]
general$cover <- speciesDominance$coverMean[match(general$taxa, speciesDominance$names)]

cor.test(general$impact, general$number)
cor.test(general$impact, general$cover)
