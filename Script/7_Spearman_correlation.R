setwd("D:/Research/Data/Redo")
library("ggpubr")
#load cut_abundance matrix and cut normalized gene data
gene <- read.csv("D:/Research/Data/Redo/log2cpm_clean_gene_rna.csv", header = TRUE,row.names = 1)
gene <-t(gene)
gene <-data.frame(gene)
genus <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus_rna.csv", header = TRUE,row.names = 1)

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

write.csv(df_cor,file="Correlation_SCC_origin_rna.csv",row.names = T)
write.csv(df_pvalue,file="Correlation_pvalue_origin_rna.csv",row.names = T)

################################################################################
####find the big correlate
temp_m <- c()
temp_g <- c()
for (i in 1:length(genus_name)) {
  for (j in 1:length(gene_name)) {
    if(abs(df_cor[[genus_name[i]]][j])>0.45){
      temp_m <- c(temp_m,genus_name[i])
      temp_g <- c(temp_g,gene_name[j])
    }
  }
}
#try to plot SCC dotplot
drawcorrelation <- function(num) {
  highscc <- cbind(genus[[temp_m[num]]],gene[[temp_g[num]]])
  highscc <- data.frame(highscc)
  
  
  temp_genus_name <- strsplit(temp_m[num],split = "\\.")
  temp_genus_name <- temp_genus_name[[1]][length(temp_genus_name[[1]])]
  
  ggscatter(highscc, x = "X1", y = "X2", 
            add = "loess", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "spearman",
            xlab = temp_genus_name, ylab = temp_g[num])+
    labs(title="Spearman correlation")
}

pdf("SCC_dotplot.pdf") 
for (i in 1:18) {
  plot(drawcorrelation(i))
}
dev.off()



########################################################################
#Immune infiltration
immune <- read.csv("D:/Research/Data/Redo/estimation_matrix.csv", header = TRUE,row.names = 1)
immune <- data.frame(t(immune))
genus <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus.csv", header = TRUE,row.names = 1)
#Sort sample name
genus <- genus[ order(row.names(genus)), ]
immune <- immune[ order(row.names(immune)), ]
mean_check <- apply(immune, 2, mean)
hist(mean_check)
sum(mean_check > 0)
immune <- immune[, mean_check > 0]
#create a blank dataframe
immune_name <- colnames(immune)
genus_name <- colnames(genus)
df_cor <- data.frame(immune_name)
df_pvalue <- data.frame(immune_name)
fill <- rep(0L, times = length(df_cor$immune_name))


for (i in genus_name) {
  df_cor$x <- fill
  df_pvalue$x<- fill
  names(df_cor)[names(df_cor) == "x"] <- i
  names(df_pvalue)[names(df_pvalue) == "x"] <- i
}
rownames(df_cor)<- df_cor$immune_name
df_cor = df_cor[,!(names(df_cor) %in% c("immune_name"))]
rownames(df_pvalue)<- df_pvalue$immune_name
df_pvalue = df_pvalue[,!(names(df_pvalue) %in% c("immune_name"))]
#calculate the correlation and fill the matrix
for (i in 1:length(genus_name)) {
  for (j in 1:length(immune_name)) {
    m<-genus_name[i]
    g<-immune_name[j]
    test<-cor.test(genus[[m]], immune[[g]], method=c("spearman"),exact=F)
    rho <- test$estimate[["rho"]]
    pvalue <- test$p.value
    df_cor[which(rownames(df_cor)==g),which(names(df_cor)==m)] <- rho
    df_pvalue[which(rownames(df_pvalue)==g),which(names(df_pvalue)==m)] <- pvalue
  }
}

write.csv(df_cor,file="Immune_Correlation_SCC_origin.csv",row.names = T)
write.csv(df_pvalue,file="Immune_Correlation_pvalue_origin.csv",row.names = T)
