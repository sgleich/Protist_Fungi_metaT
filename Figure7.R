### Protist + Fungi MetaT - Protistan metabolism over time ###
### This script will allow us to visualize the relative expression specific KO terms for ciliates and discobids over time ###
### Updated: June 11, 2024 ###

# Load libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL/FINAL")

# Load data
dfFin <- read.csv("CPM_AllEuk.csv",header=TRUE,row.names=1)


getKOs <- function(df,tax){
  # Subset dataframe for taxon of interest
  dfTax <- subset(df,grepl(paste(tax),df$Taxonomy))
  dfTax$Name <- NULL
  dfTax$Taxonomy <- NULL
  dfTax$CAZy <- NULL
  dfTax <- subset(dfTax,KEGG!="")
  
  # Sum CPM values for each KO
  dfTax <- dfTax %>% group_by(KEGG) %>% summarize_all(sum) %>% as.data.frame()
 
  # Melt dataframe
  dfMelt <- melt(dfTax,id.vars="KEGG")
  
  # Separate sample name and experiment number into different columns; Rename water column samples 
  dfNames <- colsplit(dfMelt$variable,"_",c("Sample","Experiment"))
  dfNames$Sample <- ifelse(grepl("Water",dfNames$Sample),"Water Column",dfNames$Sample)
  dfMelt <- cbind(dfMelt,dfNames)
  dfMelt$variable <- NULL
  
  # Find the KO terms that are not found in any samples from an entire experiment
  zeroKEGG <- dfMelt %>% group_by(Experiment,KEGG) %>% summarize(s=sum(value)) %>% filter(s==0) %>% as.data.frame()
  zeroKEGG <- unique(zeroKEGG$KEGG) 
  
  # Subset dataframe so that all KO terms are found in at least one sample per experiment
  `%ni%` <- Negate(`%in%`)
  dfMelt <- subset(dfMelt, KEGG %ni% zeroKEGG)
  
  # Calculate mean and sd for each KO for each experiment; will be used in the z-score standardization
  dfMean <- dfMelt %>% group_by(KEGG,Experiment) %>% summarize(mz=mean(value),sdz=sd(value)) %>% as.data.frame()
  dfMelt <- left_join(dfMelt,dfMean)
  
  # Calculate z-score
  dfMelt$z <- (dfMelt$value - dfMelt$mz)/dfMelt$sdz
  # range(dfMelt$z)
  
  # Calculate mean z-score for each KEGG, Sample, Exp (i.e., average water column duplicates)
  dfMelt <- dfMelt[c(1,3:4,7)]
  dfMeltAvg <- dfMelt %>% group_by(Experiment,Sample,KEGG) %>% summarize(mz=mean(z)) %>% as.data.frame()
  
  # List of all KO terms in dataframe to loop through
  keggTerms <- unique(dfMeltAvg$KEGG)
  
  # Establish empty dataframe
  out <- NULL
  
  # Loop through dataframe. For each KO term and experiment, find the Sample with the highest z-score.
  for (kegg in 1:length(keggTerms)){
    s <- subset(dfMeltAvg,KEGG==keggTerms[kegg])
    maxVals <- s %>% group_by(Experiment) %>% filter(mz==max(mz))
    if(length(unique(maxVals$Sample))==1){
      out <- rbind(out,maxVals)
    }
  }
  
  # Find all KO terms with max z-scores at specific time points (dependent on the temporal profile of taxa)
  if(tax=="Discob"){
    out <- subset(out,Sample=="RotT3") # Max discobid signal at T3
  }
  
  if(tax=="Cilioph"){
    out <- subset(out,Sample=="RotT6") # Max ciliate signal at T6
  }
  
  # Subset full dataset for only KO terms in out df
  subOut <- subset(dfMeltAvg,KEGG %in% out$KEGG)
  return(subOut)
}

# Run get KOs function
koOutCiliate <- getKOs(dfFin,"Cilioph")
koOutDiscob <- getKOs(dfFin,"Discob")

# Find KO terms in both ciliate and discobid dataframe
bothKOCil <- subset(koOutCiliate,KEGG %in% koOutDiscob$KEGG)
bothKODis <- subset(koOutDiscob,KEGG %in% bothKOCil$KEGG)
#unique(bothKOCil$KEGG)

# Factor levels
bothKOCil$Sample <- factor(bothKOCil$Sample,levels=c("Water Column","RotT0","RotT3","RotT6"))
bothKODis$Sample <- factor(bothKODis$Sample,levels=c("Water Column","RotT0","RotT3","RotT6"))

# Turn factors into numeric for plotting
bothKOCil$Num <- as.numeric(as.factor(bothKOCil$Sample))
bothKODis$Num <- as.numeric(as.factor(bothKODis$Sample))

# Plot function
plotFxn <- function(df,taxName){
  df$KEGGName <- ifelse(df$KEGG=="K02866","large subunit ribosomal protein L10e",NA)
  df$KEGGName <- ifelse(df$KEGG=="K02941","large subunit ribosomal protein LP0",df$KEGGName)
  df$KEGGName <- ifelse(df$KEGG=="K07375","tubulin beta",df$KEGGName)
  df$KEGGName <- ifelse(df$KEGG=="K10355","actin, other eukaryote",df$KEGGName)
  
  p <- df %>% ggplot(aes(x=Num,y=mz))+geom_point(aes(color=Experiment),size=3)+facet_wrap(~KEGGName,scales='free')+stat_smooth(method = lm, formula = y ~ poly(x, 2), se = TRUE,color="black")+scale_x_continuous(labels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))+theme_bw(base_size=16)+theme(axis.text.x = element_text(angle = 45, hjust=1), plot.title = element_text(hjust = 0.5))+xlab("Sample")+ylab(" ")+scale_color_manual(values=c("indianred","darkgreen","dodgerblue"),labels=c("Experiment #1 - 2021","Experiment #2 - 2021","Experiment #3 - 2022"))+ggtitle(paste(taxName))
  return(p)}

# Run plot function
plotOutCil <- plotFxn(bothKOCil,"Ciliate")
plotOutDis <- plotFxn(bothKODis,"Discobid")

plotOutDis + plotOutCil + plot_layout(guides = "collect",nrow=1)
ggsave("../../Figure7.pdf",width=20,height=8)
