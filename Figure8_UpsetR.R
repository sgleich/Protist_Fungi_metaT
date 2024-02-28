### Protist + Fungi MetaT - Figure 8 - UpsetR Plot ###
### This script will compare the presence/absence of CAZy group expression between fungi and all other eukaryotes ###
### Updated: February 28, 2024 ###

# Libraries
library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggupset)
library(stringi)
library(patchwork)

# Load in data
dfFin <- read.csv("CPM_all_CAZy.csv",header=TRUE,row.names=1)

# Distinguish fungi from all other eukaryotes
dfFin$Tax <- ifelse(grepl("Fungi",dfFin$Taxonomy),"Fungi","Other Eukaryote")

# Upset plot set up
dfFin$Name <- NULL
dfFin$Taxonomy <- NULL

dfMelt <- melt(dfFin,id.vars=c("CAZy","Tax"))
namez <- colsplit(dfMelt$variable,"_",c("Sample","Exp"))
namez$Sample <- ifelse(grepl("Water",namez$Sample),"Water Column",namez$Sample)
dfMelt <- cbind(dfMelt,namez)
dfMelt$variable <- NULL
dfMelt$Exp <- NULL
dfMelt <- dfMelt %>% group_by(Sample,CAZy,Tax) %>% summarize(s=sum(value)) %>% as.data.frame()
dfMelt <- subset(dfMelt,CAZy!="")

dfMelt <- dfMelt %>% pivot_wider(id_cols=c("CAZy","Sample"),names_from="Tax",values_from="s",values_fill = 0)%>% as.data.frame()

rownames(dfMelt) <- paste(dfMelt$CAZy,dfMelt$Sample,sep="_")
dfMelt$CAZy <- NULL
dfMelt$Sample <- NULL
dfMelt <- subset(dfMelt,rowSums(dfMelt)>0)


dfMeltBin <- ifelse(dfMelt > 0, 1,0)
dfMeltBin <- as.data.frame(dfMeltBin)
dfMeltBin$Intersect <- apply(dfMeltBin > 0, 1, function(x){toString(names(dfMeltBin)[x])})
dfMeltBin$Intersect <- stri_replace_all_fixed(dfMeltBin$Intersect, " ", "")
dfMeltBin$Intersect <- as.list(strsplit(dfMeltBin$Intersect, ","))

keep_df <- dfMeltBin
#keep <- list(c("Fungi"),c("Amoebozoa", "Bicoecea", "Heterolobosea","Labyrinthulomycetes"))
#keep_df <- subset(dfMeltBin,Intersect %in% keep)

# CAZy classes
keep_df$C <- ifelse(grepl("AA",rownames(keep_df)),"Auxiliary Activities",NA)
keep_df$C <- ifelse(grepl("CBM",rownames(keep_df)),"Carbohydrate-Binding Modules",keep_df$C)
keep_df$C <- ifelse(grepl("GH",rownames(keep_df)),"Glycoside Hydrolases",keep_df$C)
keep_df$C <- ifelse(grepl("GT",rownames(keep_df)),"GlycosylTransferase",keep_df$C)
keep_df$C <- ifelse(grepl("CE",rownames(keep_df)),"Carbohydrate Esterases",keep_df$C)
keep_df$C <- ifelse(grepl("PL",rownames(keep_df)),"Polysaccharide Lyases",keep_df$C)

colrs <- c("Glycoside Hydrolases"="dodgerblue",
           "GlycosylTransferase"="indianred",
           "Carbohydrate Esterases"="darkgoldenrod3",
           "Carbohydrate-Binding Modules"="forestgreen",
           "Auxiliary Activities"="salmon","Polysaccharide Lyases"="purple")

# Plot each net trap sample individually
net6 <- keep_df%>% filter(grepl("RotT6",rownames(keep_df))) %>%ggplot(aes(x=Intersect,fill=C)) + geom_bar(stat = "count", position="stack",color="black") + scale_x_upset()+theme_classic(base_size = 16)+xlab("")+ylab("Number of Shared\nand Unique CAZymes")+theme(plot.margin = margin(10, 10, 10, 100))+ggtitle(expression("Net Trap T6"))+scale_fill_manual(name=c("CAZyme Class"),values=c(colrs))+theme(legend.position="right")+ylim(0,70)


net0+net3+net6+plot_layout(guides="collect",nrow=3)+plot_annotation(tag_levels="a")
ggsave("Figure7_Feb2024.pdf",width=9,height=13)

####

fun <- subset(keep_df,Intersect=="Fungi")
fun$name <- rownames(fun)
cols <- colsplit(rownames(fun),"_",c("CAZy","Exp"))
cols %>% group_by(CAZy) %>% tally() %>% arrange(desc(n))
