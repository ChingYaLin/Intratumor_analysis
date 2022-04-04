####This Script want to cut low expression of genus normalized data and gene expression data
setwd("D:/Research/Data")
library("edgeR")
genus_raw <- read.table("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1,sep = ",")

##Clean the contaminant column
genus_raw <- genus_raw[,-c(1285:1287)]

##Calculate the mean of normalized data
mean_normalized_data <- apply(genus_raw, 2, mean)
hist(x=mean_normalized_data,
     main="Histogram of genus_Abunance",         
     xlab="Mean of normalized data",                      
     ylab="Frequency")

expr_cutoff <- 2
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(mean_normalized_data > expr_cutoff)
cut_genus <- genus_raw[, mean_normalized_data > expr_cutoff]


##deal the gene expression data
gene_raw <- read.table("D:/Research/Data/Raw/LUAD_gene_expression_raw(138sample).csv", header = TRUE,row.names = 1,sep = ",")
##Clean the non gene raws
gene_raw <- gene_raw[-c(60484:60488),]
##Remove the normal sample
gene_raw <- gene_raw[,-c(139:172)]

cpm_log <- cpm(gene_raw, log = TRUE)
cpm_log <- data.frame(cpm_log)
mean_log2_cpm <- apply(cpm_log, 1, mean)

hist(x=mean_log2_cpm,
     main="Histogram of gene expression(log2cpm)",         
     xlab="Mean of log2cpm_gene expression",                      
     ylab="Frequency")

expr_cutoff <- 3
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(mean_log2_cpm > expr_cutoff)
cut_gene <- cpm_log[mean_log2_cpm > expr_cutoff,]
cut_gene <- data.frame(cut_gene)
##rename the column by sample
name <- names(cut_gene)
for (i in 1:length(name)) {
  temp_len <- nchar(name[i])
  temp <- substr(name[i],temp_len-5,temp_len)
  names(cut_gene)[names(cut_gene) == name[i]] <- temp
}

##write the normalized gene data and genus data into CSV
write.csv(cut_gene,file="D:/Research/Data/Redo/log2cpm_clean_gene.csv",row.names = T)
write.csv(cut_genus,file="D:/Research/Data/Redo/voom_svm_clean_genus.csv",row.names = T)

#cut_gene <- rownames(cut_gene)
#cut_genus <- names(cut_genus)
