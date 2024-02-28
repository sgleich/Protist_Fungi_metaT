### Protist + Fungi MetaT - Figure 5 - KEGG expression (non-fungal eukaryotes) ###
### This script will allow us to visualize KEGG terms that are consistently the most expressed in a given sample type (i.e., Water column, net trap T0, net trap T3, or net trap T6) across the 3 decomposition experiments ###
### Updated: February 28, 2024 ###

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL")
dfFin <- read.csv("CPM_NoFungi_KEGG.csv",header=TRUE,row.names=1)

dfFin$Name <- NULL
dfFin$Taxonomy <- NULL

dfMelt <- melt(dfFin,id.vars="KEGG")
dfMelt <- dfMelt %>% filter(KEGG!="") %>%group_by(variable,KEGG) %>% summarize(s=sum(value)) %>% as.data.frame()
colz <- colsplit(dfMelt$variable,"_",c("Sample","Exp"))
colz$Sample <- ifelse(grepl("WaterColumn",colz$Sample),"WaterColumn",colz$Sample)
dfMelt <- cbind(dfMelt,colz)
dfMelt$variable <- NULL
dfMelt <- dfMelt %>% group_by(KEGG,Sample,Exp) %>% summarize(m=mean(s)) %>% as.data.frame()

subz <- subset(dfMelt, m==0)
`%ni%` <- Negate(`%in%`)
dfMelt <- subset(dfMelt,KEGG %ni% subz$KEGG)

dfMean <- dfMelt %>% group_by(KEGG,Exp) %>% summarize(mz=mean(m),sdz=sd(m)) %>% as.data.frame()
dfMelt <- left_join(dfMelt,dfMean)
dfMelt$z <- (dfMelt$m - dfMelt$mz)/dfMelt$sdz
range(dfMelt$z)

keggz <- unique(dfMelt$KEGG)
out <- NULL
for (i in 1:length(keggz)){
  sub <- subset(dfMelt,KEGG==keggz[i])
  sub2 <- sub %>% group_by(Exp) %>% top_n(1,z)
  if (sub2$Sample[1]==sub2$Sample[2] & sub2$Sample[2]==sub2$Sample[3]){
    yes <- data.frame(kegg=sub2$KEGG[1],sample=sub2$Sample[1])
    out <- rbind(out,yes)
  }
  
}

colnames(out) <- c("KEGG","Top")
dfMeltAll <- left_join(dfMelt,out)
dfMeltAll <- subset(dfMeltAll,!is.na(Top))


dfMeltAll$Sample <- ifelse(dfMeltAll$Sample=="RotT0","Net Trap T0",dfMeltAll$Sample)
dfMeltAll$Sample <- ifelse(dfMeltAll$Sample=="RotT3","Net Trap T3",dfMeltAll$Sample)
dfMeltAll$Sample <- ifelse(dfMeltAll$Sample=="RotT6","Net Trap T6",dfMeltAll$Sample)
dfMeltAll$Sample <- ifelse(dfMeltAll$Sample=="WaterColumn","Water Column",dfMeltAll$Sample)

dfMeltAll$Exp <- ifelse(dfMeltAll$Exp=="Exp1","Experiment #1 - 2021",dfMeltAll$Exp)
dfMeltAll$Exp <- ifelse(dfMeltAll$Exp=="Exp2","Experiment #2 - 2021",dfMeltAll$Exp)
dfMeltAll$Exp <- ifelse(dfMeltAll$Exp=="Exp3","Experiment #3 - 2022",dfMeltAll$Exp)

dfMeltAll$Top <- factor(dfMeltAll$Top,levels=c("WaterColumn","RotT0","RotT3","RotT6"))
dfMeltAll$Sample <- factor(dfMeltAll$Sample,levels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))
dfMeltAll$Top <- as.numeric(dfMeltAll$Top)
ggplot(dfMeltAll,aes(x=Sample,y=reorder(KEGG,Top),fill=z))+geom_tile(color="black")+facet_wrap(~Exp)+theme_classic()+scale_fill_gradient2(name="Z-score \nStandardized\nExpression",low="red",mid="white",high="blue",midpoint=0)+theme(axis.text.x = element_text(angle = 45, hjust =1))+ylab("KEGG KO")

ggsave("Figure5_Feb2024.pdf",width=7,height=12)
