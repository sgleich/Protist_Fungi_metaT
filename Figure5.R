### Protist + Fungi MetaT - WC vs. Net Trap T0 metabolism ###
### This script will allow us to visualize the relative expression of various kegg pathways in the water column and net trap T0 samples ###
### Updated: June 11, 2024 ###

# Load libraries
library(KEGGREST)
library(tidyverse)
library(stringr)
library(ggplot2)
library(reshape2)

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL/FINAL")

# Load data
dfFin <- read.csv("CPM_AllEuk.csv",header=TRUE,row.names=1)

# Load KEGG pathway info
motorp <- KEGGREST::keggLink("ko","ko04814")
motorp <- data.frame(KEGG=c(motorp),Type=c("Motor proteins (ko04814)"))
cfix <- KEGGREST::keggLink("ko","ko00710")
cfix <- data.frame(KEGG=c(cfix),Type=c("C fixation (ko00710)"))
photo <- KEGGREST::keggLink("ko","ko00195")
photo <- data.frame(KEGG=c(photo),Type=c("Photosynthesis (ko00195)"))
photoant <- KEGGREST::keggLink("ko","ko00196")
photoant <- data.frame(KEGG=c(photoant),Type=c("Photosynthesis - antenna proteins (ko00196)"))
phago<- KEGGREST::keggLink("ko","ko04145")
phago <- data.frame(KEGG=c(phago),Type=c("Phagosome (ko04145)"))
lyso <- KEGGREST::keggLink("ko","ko04142")
lyso <- data.frame(KEGG=c(lyso),Type=c("Lysosome (ko04142)"))
cmet <- KEGGREST::keggLink("ko","ko01200")
cmet <- data.frame(KEGG=c(cmet),Type=c("C Metabolism (ko01200)"))
ribo <- KEGGREST::keggLink("ko","ko03010")
ribo <- data.frame(KEGG=c(ribo),Type=c("Ribosome (ko03010)"))
protdig <- KEGGREST::keggLink("ko","ko04974")
protdig <- data.frame(KEGG=c(protdig),Type=c("Protein digestion + absorption (ko04974)"))
endo <- KEGGREST::keggLink("ko","ko04144")
endo <- data.frame(KEGG=c(endo),Type=c("Endocytosis (ko04144)"))
all <- rbind(motorp,cfix,photo,photoant,phago,lyso,cmet,ribo,protdig,endo)
all$KEGG <- str_remove_all(all$KEGG,"ko:")

# Add CPM for each KEGG term and sample
dfFin <- dfFin[c(3,5:19)]
dfFin <- dfFin %>% group_by(KEGG) %>% summarize_all(sum) %>% as.data.frame()
dfFin <- subset(dfFin,KEGG!="")

# Melt dataframe and join KEGG CPM data to KEGG pathway data
dfMelt <- melt(dfFin,id.vars=c("KEGG"))
dfMelt <- left_join(dfMelt,all)
dfMelt <- subset(dfMelt,!is.na(Type))

# Sum CPM for each pathway
dfMelt <- dfMelt %>% group_by(Type,variable)%>% summarize(s=sum(value))%>% as.data.frame()
cols <- colsplit(dfMelt$variable,"_",c("Sample","Experiment"))
cols$Sample <- ifelse(grepl("Water",cols$Sample),"Water",cols$Sample)
dfMelt <- cbind(dfMelt,cols)
dfMelt$variable <- NULL
# Average water column CPM
dfMelt <- dfMelt %>% group_by(Sample,Experiment,Type) %>% summarize(m=mean(s)) %>% as.data.frame()

# We only care about water column and net trap T0 samples
dfMelt <- subset(dfMelt,Sample=="RotT0"|Sample=="Water")

# Calculate CPM total associated with each KEGG pathway
sums <- dfMelt %>% group_by(Experiment,Type) %>% summarize(sumz=sum(m))
dfMelt <- left_join(dfMelt,sums)
dfMelt$frac <- dfMelt$m/dfMelt$sumz
#dfMelt$frac <- scales::percent(dfMelt$frac)

# Rename samples and experiments
dfMelt$Sample <- ifelse(dfMelt$Sample=="RotT0","Net Trap T0","Water Column")
dfMelt$Sample <- factor(dfMelt$Sample,levels=c("Water Column","Net Trap T0"))
dfMelt$Experiment <- ifelse(dfMelt$Experiment=="Exp1","Experiment #1 - 2021",dfMelt$Experiment)
dfMelt$Experiment <- ifelse(dfMelt$Experiment=="Exp2","Experiment #2 - 2021",dfMelt$Experiment)
dfMelt$Experiment <- ifelse(dfMelt$Experiment=="Exp3","Experiment #3 - 2022",dfMelt$Experiment)

ggplot(dfMelt,aes(x=Sample,Type,fill=frac))+geom_tile(color="black")+facet_wrap(~Experiment)+scale_fill_gradient2(name="Percent Total CPM\n(Water Column + Net Trap T0)",low="red",mid="white",high="blue",midpoint=0.5,labels = scales::percent(0.25*0:4))+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust=1))+ylab("KEGG Pathway")
ggsave("../../Figure6.pdf",width=10,height=6)
