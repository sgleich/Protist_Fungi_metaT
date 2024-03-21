# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL")

dfFin <- read.csv("CPM_all_KEGG.csv",header=TRUE,row.names=1)
dfFin <- dfFin[c(3:18)]


dfFin <- dfFin %>% filter(KEGG!="") %>% group_by(KEGG) %>% summarize_all(sum) %>% as.data.frame()

dfMelt <- melt(dfFin,id.vars="KEGG")
cols <- colsplit(dfMelt$variable,"_",c("Sample","Exp"))
dfMelt <- cbind(dfMelt,cols)
dfMelt$variable <- NULL
dfMelt$Sample <- ifelse(grepl("Water",dfMelt$Sample),"Water Column",dfMelt$Sample)
dfMelt <- dfMelt %>% group_by(Sample,Exp,KEGG) %>% summarize(m=mean(value)) %>% as.data.frame()

dfMelt$Sample <- paste(dfMelt$Sample,dfMelt$Exp,sep="_")
dfMelt$Exp <- NULL

dfWide <- dfMelt %>% pivot_wider(names_from="Sample",values_from ="m", id_cols="KEGG",values_fill = 0) %>% as.data.frame()

rownames(dfWide) <- dfWide$KEGG
dfWide$KEGG <- NULL
dfWide <- as.data.frame(t(dfWide))

dfBray <- vegdist(dfWide,method="bray")
pcoaOut <- pcoa(dfBray)

t <- data.frame(pcoaOut$vector)
t$M <- rownames(t)
cols <- colsplit(t$M,"_",c("Sample","Experiment"))
cols$Sample <- ifelse(cols$Sample=="RotT0","Net Trap T0",cols$Sample)
cols$Sample <- ifelse(cols$Sample=="RotT3","Net Trap T3",cols$Sample)
cols$Sample <- ifelse(cols$Sample=="RotT6","Net Trap T6",cols$Sample)
cols$Experiment <- ifelse(cols$Experiment=="Exp1","Experiment #1 - 2021",cols$Experiment)
cols$Experiment <- ifelse(cols$Experiment=="Exp2","Experiment #2 - 2021",cols$Experiment)
cols$Experiment <- ifelse(cols$Experiment=="Exp3","Experiment #3 - 2022",cols$Experiment)

t <- cbind(t,cols)
t$Sample <- factor(t$Sample,levels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))

pcoaOut$values

ggplot(t,aes(x=Axis.1,y=Axis.2,fill=Sample,shape=Experiment))+geom_point(size=4)+scale_fill_manual(values=c("grey","dodgerblue","maroon","darkgreen"))+scale_shape_manual(values=c(23,22,24))+ guides(fill = guide_legend(override.aes=list(shape=21)))+theme_classic()+geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept=0,linetype="dashed")+xlab("37.4%")+ylab("17.8%")
ggsave("../Figure4_March2024.pdf")
