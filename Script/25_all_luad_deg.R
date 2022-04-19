setwd("D:/Research/TCGA-LUAD all case")
library("edgeR")
library("pheatmap")
library('ggplot2')
data_raw <- read.table("D:/Research/TCGA-LUAD all case/LUAD_DEG_raw.csv", header = TRUE,row.names = 1,sep = ",")
## veiw data
dim(data_raw)
head(data_raw)
tail(data_raw)

##boxplot
data_log = log2(data_raw+1)
boxplot(data_log, cex = 1)

### remove low-expression gene
cpm_log <- cpm(data_raw, log = TRUE)
mean_log2_cpm <- apply(cpm_log, 1, mean)
hist(mean_log2_cpm)
expr_cutoff <- 0
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(mean_log2_cpm > expr_cutoff)
data_clean <- data_raw[mean_log2_cpm > expr_cutoff, ]

##check the boxplot again
data_log2 = log2(data_clean+1)
boxplot(data_log2, cex = 1,main = 'Gene expression')
cpm_log <- cpm(data_clean, log = TRUE)
mean2_data <- apply(cpm_log, 1, mean)
hist(mean2_data)

##heatmap
heatmap(cor(cpm_log),main="correlation between samples")

##pca
pca <- prcomp(t(cpm_log), scale. = TRUE)
pca_data <- pca$x
pca_data <- data.frame(pca_data)
group <- rownames(pca_data)
group <- substr(group,1,1)
pca_data$group <- group

eigs <- pca$sdev^2
pc.1 <- eigs[1] * 100 / sum(eigs)
pc.2 <- eigs[2] * 100 / sum(eigs)

pc.1 <- format(pc.1, digits = 2, nsmall = 2)
pc.2 <- format(pc.2, digits = 2, nsmall = 2)


ggplot(pca_data,aes(x=PC1,y=PC2))+
  geom_point(aes(color = as.factor(group)))+
  labs(title = "PCA",x=paste("PC1 (", pc.1, "%)", sep = ""),y=paste("PC2 (", pc.2, "%)", sep = ""))+
  scale_color_discrete(name = "Group", labels = c("Normal", "Tumor"))

plot(pca$x[, 1], pca$x[, 2], pch = 21, xlab = paste("PC1 (", pc.1, "%)", sep = ""), 
     ylab = paste("PC2 (", pc.2, "%)", sep = ""),main = "PCA")
text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))


##make group(2 group:ad,n)
group <- substr(colnames(data_clean), 1, 1)
y <- DGEList(counts = data_clean, group = group)
y <- calcNormFactors(y) #TMM normalized
y$samples
y <- estimateDisp(y) 
sqrt(y$common.dispersion) # biological coefficient of variation
plotBCV(y)

#find deg
et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(data_clean), sort.by = "none")
head(results_edgeR$table)
write.csv(results_edgeR$table,file="DEG_result.csv",row.names = T)
#{r count-de-edgeR}
#How many genes are differentially expressed at an FDR of 10%?
#MA-plot
sum(results_edgeR$table$FDR < 1e-10)
sum(results_edgeR$table$FDR < 1e-5)
plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < 1e-10],main = 'MA-plot')
abline(h = 0, col = "blue")

#volcano plot
volcanoData <- cbind(results_edgeR$table$logFC, -log10(results_edgeR$table$FDR))
test <- data.frame(volcanoData)
test$color <- "no"
test$color[test$X1>2&test$X2>10] <- "upregulate"
test$color[test$X1<(-2)&test$X2>10] <- "downregulate"
colnames(test) <- c("logFC", "logFDR","color")
head(volcanoData)
ggplot(data=test, aes(x=logFC, y=logFDR, col=color)) +
  geom_point(size =1) + 
  theme_minimal() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=10, col="red")+
  labs(title = "Volcano plot")


#select gene and plot heatmap
selY <- data_clean[rownames(results_edgeR$table)[results_edgeR$table$PValue< 0.05 & 
                                                   abs(results_edgeR$table$logFC)>2],]

cpm_log2_sel <- cpm(selY, log = TRUE)
cpm_log2_sel <- data.frame(cpm_log2_sel)


sely_log <- log10(selY+1)
head(sely_log)
selected_list <- list(rownames(selY))

heatmap(t(t(sely_log)),scale = 'row',main='heatmap')

write.csv(selY,file="select_gene_raw.csv",row.names = T)
write.csv(cpm_log2_sel,file="select_gene_log2cpm.csv",row.names = T)
write.csv(selected_list[[1]],file = "DEG_namelist.csv",row.names = F)

####Check overlap with origin deg analysis####
all_deg <- selected_list[[1]]
my_deg <- read.csv("D:/Research/Data/DEG_namelist.csv",sep = "\t",header = T)
my_deg <- my_deg[order(my_deg$x),]


library(ggvenn)
x<- list(All_sample = all_deg,
         Data_sample = my_deg)

ggvenn(
  x, 
  fill_color = c("#FFCC99", "#99CCFF"),
  stroke_color = "#C0C0C0",
  stroke_size = 1, set_name_size = 8,
  label_sep = ",",
  text_size = 5
)