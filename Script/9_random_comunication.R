setwd("D:/Research/Data/Redo")
gene <- read.csv("D:/Research/Data/Redo/log2cpm_clean_gene_rna.csv", header = TRUE,row.names = 1)
gene <-t(gene)
gene <-data.frame(gene)
genus <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus_rna.csv", header = TRUE,row.names = 1)

#gene <- gene[ order(row.names(gene)), ]
#mean_check <- apply(gene, 2, mean)
#hist(mean_check)
#sum(mean_check > 0)
#gene <- gene[, mean_check > 0]
#random communication(reorder sample value)
for (i in 1:length(colnames(gene))) {
  temp <- colnames(gene)[i]
  gene[[temp]] <-sample(gene[[temp]])
}

#sort the row name
genus <- genus[ order(row.names(genus)), ]
gene <- gene[ order(row.names(gene)), ]

#create a blank dataframe
gene_name <- colnames(gene)
genus_name <- colnames(genus)
df_cor <- data.frame(gene_name)
df_pvalue <- data.frame(gene_name)
fill <- rep(0L, times = length(df_cor$gene_name))

for (i in genus_name) {
  df_cor$x <- fill
  df_pvalue$x<- fill
  names(df_cor)[names(df_cor) == "x"] <- i
  names(df_pvalue)[names(df_pvalue) == "x"] <- i
}
rownames(df_cor)<- df_cor$gene_name
df_cor = df_cor[,!(names(df_cor) %in% c("gene_name"))]
rownames(df_pvalue)<- df_pvalue$gene_name
df_pvalue = df_pvalue[,!(names(df_pvalue) %in% c("gene_name"))]

#calculate the correlation and fill the matrix
for (i in 1:length(genus_name)) {
  for (j in 1:length(gene_name)) {
    m<-genus_name[i]
    g<-gene_name[j]
    test<-cor.test(genus[[m]], gene[[g]], method=c("spearman"),exact=F)
    rho <- test$estimate[["rho"]]
    pvalue <- test$p.value
    df_cor[which(rownames(df_cor)==g),which(names(df_cor)==m)] <- rho
    df_pvalue[which(rownames(df_pvalue)==g),which(names(df_pvalue)==m)] <- pvalue
  }
}

write.csv(df_cor,file="Correlation_SCC_origin_rna_random.csv",row.names = T)
write.csv(df_pvalue,file="Correlation_pvalue_origin_rna_random.csv",row.names = T)



