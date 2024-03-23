### Protist + Fungi MetaT - Figure 5 - KEGG expression (non-fungal eukaryotes) ###
### This script will allow us to visualize KEGG terms that are consistently the most expressed in a given sample type (i.e., Water column, net trap T0, net trap T3, or net trap T6) across the 3 decomposition experiments ###
### Updated: February 28, 2024 ###

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL")
dfFin <- read.csv("CPM_All_KEGG.csv",header=TRUE,row.names=1)

dfFin$Name <- NULL
dfFin$Taxonomy <- NULL

dfMelt <- melt(dfFin,id.vars="KEGG")
dfMelt <- dfMelt %>% filter(KEGG!="") %>%group_by(variable,KEGG) %>% summarize(s=sum(value)) %>% as.data.frame()
colz <- colsplit(dfMelt$variable,"_",c("Sample","Exp"))
colz$Sample <- ifelse(grepl("WaterColumn",colz$Sample),"WaterColumn",colz$Sample)
dfMelt <- cbind(dfMelt,colz)
dfMelt$variable <- NULL
dfMelt <- dfMelt %>% group_by(KEGG,Sample,Exp) %>% summarize(m=mean(s)) %>% as.data.frame()

dfMelt <- subset(dfMelt,Sample=="RotT0"|Sample=="WaterColumn")


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

dfMeltAll <- dfMeltAll[c(1,2,3,4)]

dfWide <- dfMeltAll %>% pivot_wider(id_cols=c("KEGG","Exp"),names_from="Sample",values_from = "m" )

e1 <- subset(dfWide,Exp=="Experiment #1 - 2021")
e2 <- subset(dfWide,Exp=="Experiment #2 - 2021")
e3 <- subset(dfWide,Exp=="Experiment #3 - 2022")

e1$col <- ifelse(e1$`Net Trap T0`>e1$`Water Column`,"dodgerblue","grey")
e2$col <- ifelse(e2$`Net Trap T0`>e2$`Water Column`,"dodgerblue","grey")
e3$col <- ifelse(e3$`Net Trap T0`>e3$`Water Column`,"dodgerblue","grey")

e4 <- e3
e4$`Net Trap T0` <- 0
e4$`Water Column` <- 0

exp4<- ggplot(e4,aes(x=`Water Column`,y=`Net Trap T0`,fill=col))+theme_bw(base_size=12)+geom_abline (slope=1, linetype = "dashed", color="black")+scale_fill_manual(name="Sample",breaks=c("grey","dodgerblue"),values=c("grey","dodgerblue"),labels=c("Water Column","Net Trap T0"))+xlab("Water Column Mean CPM")+ylab("Net Trap T0 CPM")+ggtitle("Top KO Terms")+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+scale_x_continuous(labels = function(x) format(x, scientific = TRUE))#+geom_point(shape=21,size=2)

exp1+exp2+exp3+exp4+plot_layout(guides = "collect",nrow=2)
ggsave("../WC_RotT0.pdf",width=10,height=8)

e1a <- e1 %>% filter(col=="dodgerblue") %>% arrange(desc(`Net Trap T0`))
e1a <- e1a[1:20,]
e2a <- e2 %>% filter(col=="dodgerblue") %>% arrange(desc(`Net Trap T0`))
e2a <- e2a[1:20,]
e3a <- e3 %>% filter(col=="dodgerblue") %>% arrange(desc(`Net Trap T0`))
e3a <- e3a[1:20,]
keepa <- subset(e1a, KEGG %in% e2a$KEGG & KEGG %in% e3a$KEGG)

e1b <- e1 %>% filter(col=="grey") %>% arrange(desc(`Water Column`))
e1b <- e1b[1:20,]
e2b <- e2 %>% filter(col=="grey") %>% arrange(desc(`Water Column`))
e2b <- e2b[1:20,]
e3b <- e3 %>% filter(col=="grey") %>% arrange(desc(`Water Column`))
e3b <- e3b[1:20,]

keepb <- subset(e1b, KEGG %in% e2b$KEGG & KEGG %in% e3b$KEGG)

new <- c(unique(keepa$KEGG),unique(keepb$KEGG))

dfNew %>% filter(Exp=="Experiment #1 - 2021")%>% ggplot(aes(x=Sample,y=KEGG,fill=m))+geom_tile()
