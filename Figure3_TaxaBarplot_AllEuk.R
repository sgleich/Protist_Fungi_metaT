### Protist + Fungi MetaT - Figure 3 - Taxa barplot whole eukaryotic community ###
### This script will allow us to visualize the taxonomic breakdown to the microeukaryote community across the 3 decomposition experiments ###
### Updated: February 26, 2024 ###

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL")

dfFin <- read.csv(file.choose(),header=TRUE,row.names=1)

dfTax <- dfFin[c(2,4:18)]

taxz <- colsplit(dfTax$Taxonomy,";",c("d","k","p","c","o","f","g","s"))
taxz$fin <- ifelse(taxz$c==" Dinophyceae","Dinoflagellate",NA)
taxz$fin <- ifelse(taxz$p==" Choanoflagellida","Choanoflagellate",taxz$fin)
taxz$fin <- ifelse(grepl("MAST",taxz$p),"MAST",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Foraminifera","Foraminifera",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Haptophyta","Haptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Chlorophyta","Chlorophyte",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Cryptophyta","Cryptophyte",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Fungi","Fungi",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Cercozoa","Cercozoa",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Ciliophora","Ciliate",taxz$fin)
taxz$fin <- ifelse(taxz$p==" Discoba","Discoba",taxz$fin)
taxz$fin <- ifelse(taxz$k==" Amoebozoa","Amoebozoa",taxz$fin)
taxz$fin <- ifelse(taxz$c==" Bacillariophyta","Diatom",taxz$fin)
taxz$fin <- ifelse(taxz$c==" Pelagophyceae","Pelagophyte",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin) & taxz$k==" Stramenopiles","Other Stramenopiles",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin) & taxz$k==" Alveolata","Other Alveolate",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin) & taxz$k==" Archaeplastida","Other Archaeplastida",taxz$fin)
taxz$fin <- ifelse(is.na(taxz$fin),"Other Eukaryote",taxz$fin)
unique(taxz$fin)

dfTax$tax <- taxz$fin
dfTax$Taxonomy <- NULL
dfTaxMelt <- melt(dfTax,id.vars="tax")
dfTaxMelt <- dfTaxMelt %>% group_by(variable,tax)%>% summarize(s=sum(value))

varSplit <- colsplit(dfTaxMelt$variable,"_",c("Sample","Exp"))
dfTaxMelt <- cbind(dfTaxMelt,varSplit)
colrs <- c("Choanoflagellate"="#D46854", "Foraminifera"="#BFA88C","Haptophyte"="#993EE8","Other Stramenopiles"="#D895DD" ,"Chlorophyte"="#72DBBD","Cryptophyte"="#E6B15C" ,"Amoebozoa"="#DBE7D4","Dinoflagellate"="#D2E04E","Discoba"="#72E16B","Diatom"="#D2608C","Other Alveolate"="#DAB8CE","Other Eukaryote"="#7D98D5","Cercozoa"="#C6E198","Fungi"="#84C3D7","Ciliate"="#8169D7","Pelagophyte"="#DB57D1","Other Archaeplastida"="dodgerblue","MAST"="khaki1")

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
ggsave("Figure3_Feb2024.pdf",width=17,height=7)
