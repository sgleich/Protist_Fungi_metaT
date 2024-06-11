### Protist + Fungi MetaT - Pie chart comparison ###
### This script will allow us to visualize the differences in the number of fungal vs. non-fungal transcripts when using different databases to assign taxonomy to metaT reads. ###
### Updated: June 11, 2024 ###

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL/FINAL")

# Load datasets IDed with ProFun and MMETSP
dfMyco <- read.csv("CPM_AllEuk.csv",header=TRUE,row.names=1)
dfMMETSP <- read.csv("CPM_AllEuk_MMETSP.csv",header=TRUE,row.names=1)

# Remove columns we don't need
dfMyco$Name <- NULL
dfMyco$KEGG <- NULL
dfMyco$CAZy <- NULL
dfMMETSP$Name <- NULL
dfMMETSP$KEGG <- NULL
dfMMETSP$CAZy <- NULL

# Labels fungi vs. other eukaryotes
dfMyco$Fin <- ifelse(grepl("Fungi",dfMyco$Taxonomy),"Fungi","Other Eukaryote")
dfMMETSP$Fin <- ifelse(grepl("Fungi",dfMMETSP$Taxonomy),"Fungi","Other Eukaryote")
dfMyco$Taxonomy <- NULL
dfMMETSP$Taxonomy <- NULL

# Melt data
dfMycoM <- melt(dfMyco,id.vars="Fin")
dfMMEM <- melt(dfMMETSP,id.vars="Fin")

# Summarize data
dfMycoM <- dfMycoM %>% group_by(Fin,variable) %>% summarize(s=sum(value)) %>% as.data.frame()
dfMMEM <- dfMMEM %>% group_by(Fin,variable) %>% summarize(s=sum(value)) %>% as.data.frame()

sumzMyco <- dfMycoM %>% group_by(variable) %>% summarize(s2=sum(s))
dfMycoM <- left_join(dfMycoM,sumzMyco)
dfMycoM$per <- dfMycoM$s/dfMycoM$s2

sumzMMEM <- dfMMEM %>% group_by(variable) %>% summarize(s2=sum(s))
dfMMEM<- left_join(dfMMEM,sumzMMEM)
dfMMEM$per <- dfMMEM$s/dfMMEM$s2

# Rename sample types
dfMycoM$Sample <- ifelse(grepl("RotT0",dfMycoM$variable),"Net Trap T0","Water Column")
dfMycoM$Sample <- ifelse(grepl("RotT3",dfMycoM$variable),"Net Trap T3",dfMycoM$Sample)
dfMycoM$Sample <- ifelse(grepl("RotT6",dfMycoM$variable),"Net Trap T6",dfMycoM$Sample)

dfMMEM$Sample <- ifelse(grepl("RotT0",dfMMEM$variable),"Net Trap T0","Water Column")
dfMMEM$Sample <- ifelse(grepl("RotT3",dfMMEM$variable),"Net Trap T3",dfMMEM$Sample)
dfMMEM$Sample <- ifelse(grepl("RotT6",dfMMEM$variable),"Net Trap T6",dfMMEM$Sample)


# Calculate average fungal vs. non-fungal transcripts across the 3 experiments
dfMycoM$variable <- NULL
dfMycoM$s <- NULL
dfMycoM$s2 <- NULL
dfMycoM <- dfMycoM %>% group_by(Fin,Sample) %>% summarize(m=mean(per)) %>% as.data.frame()

dfMMEM$variable <- NULL
dfMMEM$s <- NULL
dfMMEM$s2 <- NULL
dfMMEM <- dfMMEM %>% group_by(Fin,Sample) %>% summarize(m=mean(per)) %>% as.data.frame()

# Summarize
df <- dfMMEM %>% 
  group_by(Sample) %>% # Variable to be transformed
  mutate(perc = m/ sum(m)) %>% 
  arrange(perc) %>%
  mutate(labels = scales::percent(perc))

# Factor sample types
df$Sample <- factor(df$Sample,levels=c("Water Column","Net Trap T0","Net Trap T3","Net Trap T6"))

# Plot
mmep <- ggplot(df, aes(x = "", y = perc, fill = Fin)) +
  geom_col(color="black") + geom_label(aes(label = labels), position = position_stack(vjust = 0.5), show.legend = FALSE,color="white") +
  coord_polar(theta = "y")+facet_wrap(~Sample,nrow=1,ncol=4)+theme_void()+scale_fill_manual(name="Taxonomic Group",values=c("black","grey"))+ggtitle("MMETSP")

mmep+mycop+plot_layout(guides = "collect",nrow=2)
ggsave("PieChart.pdf")
