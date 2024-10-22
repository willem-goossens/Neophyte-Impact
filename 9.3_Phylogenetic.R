
rm(list=ls())
library(V.PhyloMaker2)
library(readr)

impact <-read_csv("impact_diaz.csv")
#impact <- impact[impact$class=="0%-1%",]
#impact <- impact[!is.na(impact$LMA),]

#impact <- impact[!duplicated(impact$taxa),]

phylo <- read_csv("phylo.csv")


tree <- phylo.maker(sp.list=phylo, tree= GBOTB.extended.TPL, nodes= nodes.info.1.TPL, scenarios = "S3")
#write.tree(tree$scenario.3, "sample.tre")
# plot(tree$scenario.3)


#library(ggtree)
#ggtree(tree$scenario.3, layout="circular", ladderize = FALSE)+ geom_nodepoint() + geom_tiplab(hjust = -.1)

library(nlme)
library(ape)
library(caper)


colnames(impact)[1] <- c("Species")
impact$Species <- (phylo$species)
impact$Species <- gsub(" ","_", impact$Species)
impact <- impact[!(impact$Species=="crop_vineyard"),]


pgls_replicated <- gls(RelDiff ~ LMA + class + SSD + LA+H+SM+ N, 
                       correlation = corBrownian(1, tree$scenario.3, form= ~Species),  # Use the Brownian motion model
                       method = "ML", 
                       data = impact, na.action = na.omit)
pgls_replicated
summary(pgls_replicated)
plot(pgls_replicated)


