setwd("D:/Research/Data/Raw")
library(TCGAutils)
data <- read.csv("Metadata_LUAD_RNA_only138.csv")

uuids <- data$case_uuid

barcode <- UUIDtoBarcode(uuids, from_type = "case_id")
write.csv(barcode,file = 'TCGA_Barcode_rna.csv',row.names = F)
