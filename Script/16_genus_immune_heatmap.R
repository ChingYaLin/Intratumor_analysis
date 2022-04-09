immune_genus_scc <- read.csv("D:/Research/Data/Redo/Immune_Correlation_SCC_origin.csv", header = TRUE,row.names = 1)
library(tidyverse)
DF <- rownames_to_column(immune_genus_scc, var = "Xval")
DF <- DF %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")

select_combination <- DF[which(abs(DF$Value)>=0.3),]
select_combination <- select_combination[order(select_combination$Source),]

genus <- select_combination$Source
genus <- genus[!duplicated(genus)]
#select 45 selected genus to SCC table
plot_table <- immune_genus_scc[,which(names(immune_genus_scc)%in%genus)]
#rename the genus name
colname_genus <- colnames(plot_table)
for (j in 1:length(colname_genus)) {
  temp <- strsplit(colname_genus[j],split = "\\.")
  names(plot_table)[names(plot_table)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
}

library("pheatmap")
plot_table <- plot_table[28:48,]
pdf(file = "D:/Research/Data/Visualization/16.GenusImmuneHeatmap.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 30)
pheatmap(plot_table,cutree_cols  = 3,main = "The correlation of single genus and immune infiltration")
dev.off()

genus_abundance <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)
genus_abundance <- genus_abundance[,which(names(genus_abundance)%in%genus)]
#rename the genus name
colname_genus <- colnames(genus_abundance)
for (j in 1:length(colname_genus)) {
  temp <- strsplit(colname_genus[j],split = "\\.")
  names(genus_abundance)[names(genus_abundance)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
}

pheatmap(genus_abundance,scale = "column",
         color = hcl.colors(400, "RdBu"),
         breaks = test)
BluYl
breaks = c(2,1.5,1,0.5,0,-0.5,-1,-1.5,-2)
test <- 200:-200
test <- test*0.01
