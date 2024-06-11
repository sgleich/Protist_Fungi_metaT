### Protist + Fungi MetaT - Taxa barplot fungal community ###
### This script will allow us to visualize the taxonomic breakdown to the fungal community across the 3 decomposition experiments ###
### Updated: June 11, 2024 ###

# Load libraries
library(ggpubr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(patchwork)
library(stringr)
library(randomcoloR)

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL/FINAL")

# Load data
dfFin <- read.csv("CPM_Fungi.csv",header=TRUE,row.names=1)

# Prep data for taxonomy
dfTax <- dfFin[c(2,5:19)]
dfTax <- subset(dfTax,Taxonomy!="Eukaryota;Opisthokonta;Fungi")
dfTax <- subset(dfTax,Taxonomy!="Eukaryota; Opisthokonta; Fungi")

# Prep fungal phyla
taxz <- colsplit(dfTax$Taxonomy,";",c("d","k","p","c","o","f","g","s"))
taxz$c <- str_remove_all(taxz$c," ")
unique(taxz$c)

dfTax$tax <- taxz$c
dfTax$Taxonomy <- NULL

# Summarize data
dfTaxMelt <- melt(dfTax,id.vars="tax")
dfTaxMelt <- dfTaxMelt %>% group_by(variable,tax)%>% summarize(s=sum(value))
varSplit <- colsplit(dfTaxMelt$variable,"_",c("Sample","Exp"))
dfTaxMelt <- cbind(dfTaxMelt,varSplit)
colrs <- distinctColorPalette(length(unique(dfTaxMelt$tax)))

# Tax plot function
taxPltFxn <- function(df,exp,title){
  subs <- subset(df, Exp==paste(exp))
  subsWater <- subset(subs,grepl("Water",subs$Sample))
  subsWater <- subsWater %>% group_by(tax) %>% summarize(s=mean(s))
  subsWater$Sample <- "Water Column"
  subs <- subset(subs,!grepl("Water",subs$Sample))
  subs$variable <- NULL
  subs$Exp <- NULL
  subs$Sample <- ifelse(subs$Sample=="RotT0",0,subs$Sample)
  subs$Sample <- ifelse(subs$Sample=="RotT3",3,subs$Sample)
  subs$Sample <- ifelse(subs$Sample=="RotT6",6,subs$Sample)
  
  
  waterPlot <- ggplot(subsWater,aes(x=Sample,y=s,fill=tax))+geom_bar(stat="identity",color="grey20",position="fill")+scale_fill_manual(values=c(colrs))+theme_classic(base_size=16)+xlab("")+ylab("Relative Transcript Abundance")+theme(axis.text.x = element_text(angle = 45, hjust =1))+theme(legend.position = "none")
  
  rotPlot <- ggplot(subs,aes(x=as.numeric(Sample),y=s,fill=tax))+geom_area(alpha = 0.6, position = 'fill')+geom_col(width = 1.5, color = 'gray20', position = 'fill')+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+theme_classic(base_size = 16)+xlab("")+ylab("Relative Transcript Abundance")+theme(axis.text.x = element_text(angle = 45, hjust =1))+scale_x_continuous(breaks=c(0,3,6),labels=c("Net Trap\nDay 0","Net Trap\nDay 3","Net Trap\nDay 6"))+theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
  
  out <- waterPlot+rotPlot+plot_layout(widths=c(5,20),guides = "collect",nrow=1)+plot_annotation(title=paste(title))
  return(out)}

pExp1 <- taxPltFxn(dfTaxMelt,"Exp1","Experiment #1\n2021")
pExp2 <- taxPltFxn(dfTaxMelt,"Exp2","Experiment #2\n2021")
pExp3 <- taxPltFxn(dfTaxMelt,"Exp3","Experiment #3\n2022")
ggarrange(pExp1,pExp2,pExp3,common.legend = TRUE,nrow=1,ncol=3,labels=c("a","b","c"),font.label = list(size=12,color="black",face="plain"),hjust=-0.1,vjust=0.5)
ggsave("../../FungalBarplot.pdf",width=17,height=7)
