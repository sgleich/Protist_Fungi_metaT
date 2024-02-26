### Protist + Fungi MetaT - Library Normalization and CPM Calculation ###
### This script will normalize and calculate the CPM of the decomposition experiment metaT data for Fungi ONLY to use in downstream KEGG KO analyses ###
### Updated: February 26, 2024 ###

# Set working directory
setwd("~/Desktop/PARAGON_DATA_FINAL")

# Load in taxonomic IDs (dfEuk), functional IDs (dfEgg), and transcript abundances (dfSalmon)
dfEuk <- read.delim(file.choose(),header=TRUE,sep="\t")
dfEgg <- read.delim(file.choose(),header=TRUE,skip=3)
dfSalmon <- read.csv(file.choose(),header=TRUE,row.names=1)

# Get taxonomic ID with the max pid.
dfEuk <- dfEuk[c(2,4,6)]
dfEuk <- dfEuk %>% arrange(desc(max_pid)) %>% distinct(transcript_name,.keep_all = TRUE) %>% as.data.frame() 
dfEuk$max_pid <- NULL

# Get the eggnog-mapper columns we'll be using for functional annotations: KEGG ko and CAZy
dfEgg <- dfEgg[c(1,9,16)]
colnames(dfEgg) <- c("transcript_name","KEGG","CAZy")

# Join EUKulele and eggnog-mapper tables and modify names of proteinIDs
dfEukEgg <- left_join(dfEuk,dfEgg)
dfEukEgg <- subset(dfEukEgg,grepl(".p1",dfEukEgg$transcript_name))
colz <- colsplit(dfEukEgg$transcript_name,"\\.p",c("name","p"))
dfEukEgg$transcript_name <- colz$name
colnames(dfSalmon)[1] <- "transcript_name"

# Join EUKulele + eggnog-mapper table and salmon (transcript counts) table.
dfAll <- left_join(dfEukEgg,dfSalmon)

# Grab all eukaryotic proteinIDs - we don't care about anything that is not eukaryotic.
dfAll <- subset(dfAll,grepl("Eukaryot",dfAll$full_classification))

# Parse KEGG ko data.
dfAll$KEGG <- ifelse(is.na(dfAll$KEGG),"",dfAll$KEGG)
dfAll$CAZy <- ifelse(is.na(dfAll$CAZy),"",dfAll$CAZy)
dfAll$CAZy <- NULL
dfAll$KEGG <- str_remove_all(dfAll$KEGG,"ko:")
dfAll$KEGG <- ifelse(is.na(dfAll$KEGG),"",dfAll$KEGG)
dfAll$KEGG <- str_split(dfAll$KEGG,",")
dfWide <- unnest(dfAll,KEGG)
colnames(dfWide)[1:3] <- c("Name","Taxonomy","KEGG")
dfWide  <- dfWide %>% distinct(.,.keep_all = TRUE) %>% as.data.frame()

# Subset for fungi
dfWide <- subset(dfWide,grepl("Fungi",dfWide$Taxonomy))

# Separate Exp #1 and Exp #2 from Exp #3 (Experiment #3 was sequenced separately)
dfExp12 <- dfWide[,c("Name","Taxonomy","KEGG","RotT0_Exp1","RotT3_Exp1","RotT6_Exp1","RotT0_Exp2","RotT3_Exp2","RotT6_Exp2","WaterColumn1_Exp1","WaterColumn2_Exp1","WaterColumn1_Exp2","WaterColumn2_Exp2")]
dfExp3 <- dfWide[,c("Name","Taxonomy","KEGG","RotT0_Exp3","RotT3_Exp3","RotT6_Exp3","WaterColumn1_Exp3","WaterColumn2_Exp3")] 

# Normalize Exp 1&2 Libraries and Calculate CPM
dgeExp12_list <- DGEList(counts=dfExp12[4:13],genes=dfExp12[1:3],group=c(rep("T0_Exp1",1),rep("T3_Exp1",1),rep("T6_Exp1",1),rep("T0_Exp2",1),rep("T3_Exp2",1),rep("T6_Exp2",1),rep("Water_Exp1",2),rep("Water_Exp2",2)))
dgeExp12_list$samples # Visualize groups

Exp12Norm <- calcNormFactors(dgeExp12_list,method = "TMM")
Exp12Cpm <- cpm(Exp12Norm, normalized.lib.sizes=TRUE, log=FALSE) 
Exp12Cpm <-as.data.frame(Exp12Cpm)  
Exp12Cpm <- data.frame(Exp12Norm$genes,Exp12Cpm)

# Normalize Exp 3 Libraries and Calculate CPM 
dgeExp3 <- DGEList(counts=dfExp3[4:8],genes=dfExp3[1:3],group=c(rep("T0_Exp3",1),rep("T3_Exp3",1),rep("T6_Exp3",1),rep("Water_Exp3",2)))
dgeExp3$samples # Visualize groups

Exp3Norm <- calcNormFactors(dgeExp3,method = "TMM")
Exp3Cpm <- cpm(Exp3Norm, normalized.lib.sizes=TRUE, log=FALSE) 
Exp3Cpm <-as.data.frame(Exp3Cpm)  
Exp3Cpm <-data.frame(Exp3Norm$genes,Exp3Cpm)

# Combine Normalized Exp 1&2 data with Normalized Exp 3 data
dfFin <- left_join(Exp12Cpm,Exp3Cpm)

# Save CPM values for downstream analyses
write.csv(dfFin,"CPM_Fungi_KEGG.csv")
