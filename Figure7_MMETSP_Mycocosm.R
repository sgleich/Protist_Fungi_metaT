### Protist + Fungi MetaT - Figure 7 - Taxa barplot fungi ###
### This script will allow us to visualize the differences in the number of fungal vs. non-fungal transcripts when using different databases to assign taxonomy to metaT reads. ###
### Updated: February 26, 2024 ###

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL")

dfMyco <- read.csv(file.choose(),header=TRUE,row.names=1)
dfMMETSP <- read.csv(file.choose(),header=TRUE,row.names=1)

dfMyco$Name <- NULL
dfMyco$KEGG <- NULL
dfMMETSP$Name <- NULL
dfMMETSP$KEGG <- NULL

dfMyco$Fin <- ifelse(grepl("Fungi",dfMyco$Taxonomy),"Fungi","Other Eukaryote")
dfMMETSP$Fin <- ifelse(grepl("Fungi",dfMMETSP$Taxonomy),"Fungi","Other Eukaryote")

dfMyco$Taxonomy <- NULL
dfMMETSP$Taxonomy <- NULL

dfMycoM <- melt(dfMyco,id.vars="Fin")
dfMMEM <- melt(dfMMETSP,id.vars="Fin")

dfMycoM <- dfMycoM %>% group_by(Fin,variable) %>% summarize(s=sum(value)) %>% as.data.frame()
dfMMEM <- dfMMEM %>% group_by(Fin,variable) %>% summarize(s=sum(value)) %>% as.data.frame()


dfMycoM$Sample <- ifelse(grepl("RotT0",dfMycoM$variable),"Net Trap T0","Water Column")
dfMycoM$Sample <- ifelse(grepl("RotT3",dfMycoM$variable),"Net Trap T3",dfMycoM$Sample)
dfMycoM$Sample <- ifelse(grepl("RotT6",dfMycoM$variable),"Net Trap T6",dfMycoM$Sample)

dfMMEM$Sample <- ifelse(grepl("RotT0",dfMMEM$variable),"Net Trap T0","Water Column")
dfMMEM$Sample <- ifelse(grepl("RotT3",dfMMEM$variable),"Net Trap T3",dfMMEM$Sample)
dfMMEM$Sample <- ifelse(grepl("RotT6",dfMMEM$variable),"Net Trap T6",dfMMEM$Sample)

dfMycoM$variable <- NULL
dfMycoM <- dfMycoM %>% group_by(Fin,Sample) %>% summarize(m=mean(s)) %>% as.data.frame()

dfMMEM$variable <- NULL
dfMMEM <- dfMMEM %>% group_by(Fin,Sample) %>% summarize(m=mean(s)) %>% as.data.frame()

df <- dfMycoM %>% 
  group_by(Sample) %>% # Variable to be transformed
  mutate(perc = m/ sum(m)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

df$Sample <- factor(df$Sample,levels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))

mycop <- ggplot(df, aes(x = "", y = perc, fill = Fin)) +
  geom_col(color="black") + geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE,color="white") +
  coord_polar(theta = "y")+facet_wrap(~Sample,nrow=1,ncol=4)+theme_void()+scale_fill_manual(name="Taxonomic Group",values=c("black","grey"))+ggtitle("MMETSP + Mycocosm")

mmep+mycop+plot_layout(guides = "collect",nrow=2)
ggsave("Figure7_Feb2024.pdf")
