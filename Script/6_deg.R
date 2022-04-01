setwd("D:/Research/Data")
library("edgeR")
library("pheatmap")
library('ggplot2')

data_raw <- read.table("D:/Research/Data/Raw/LUAD_gene_expression_raw(138sample).csv", header = TRUE,row.names = 1,sep = ",")
## veiw data
dim(data_raw)
head(data_raw)
tail(data_raw)
##remove first column because it is gene name
data_clean <- data_raw
data_log = log2(data_clean+1)
boxplot(data_log, cex = 1)

### remove outlier(remove expression 0 first)
cpm_log <- cpm(data_clean, log = TRUE)
mean_log2_cpm <- apply(cpm_log, 1, mean)
hist(mean_log2_cpm)
expr_cutoff <- 1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
data_clean <- data_clean[median_log2_cpm > expr_cutoff, ]


##check the boxplot again
data_log2 = log2(data_clean+1)
boxplot(data_log2, cex = 1,main = 'Gene expression')
cpm_log <- cpm(data_clean, log = TRUE)
mean2_data <- apply(cpm_log, 1, mean)
hist(mean2_data)

##heatmap
heatmap(cor(cpm_log),main="correlation between samples")
pheatmap(cor(cpm_log),main = "correlation between samples")

##pca
pca <- prcomp(t(cpm_log), scale. = TRUE)
pca_data <- pca$x
pca_data <- data.frame(pca_data)
group <- rownames(pca_data)
group <- substr(group,1,2)
pca_data$group <- group
ggplot(pca_data,aes(x=PC1,y=PC2))+
  geom_point(aes(color = as.factor(group)))+
  labs(title = "PCA",x="PC1(50.8125)",y="PC2(47.0107)")+
  scale_color_discrete(name = "Group", labels = c("Normal", "Pair_Tumor", "Tumor"))

plot(pca$x[, 1], pca$x[, 2], pch = 21, xlab = "PC1(50.8125)", ylab = "PC2(47.0107)",main = "PCA")
text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
#summary(pca)

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
#remove the normal column
cpm_log2_sel <- cpm_log2_sel[,-c(139:172)]
#rename the column by sample
name <- names(cpm_log2_sel)
for (i in 1:length(name)) {
  temp_len <- nchar(name[i])
  temp <- substr(name[i],temp_len-5,temp_len)
  names(cpm_log2_sel)[names(cpm_log2_sel) == name[i]] <- temp
}

sely_log <- log10(selY+1)
head(sely_log)
selected_list <- list(rownames(selY))

pheatmap(sely_log,
         scale = "row",
         main = 'Heatmap',
         color = colorRampPalette(c("red","white", "blue"))(4),
         breaks = c(2, 1,0,-1,-2))
heatmap(t(t(sely_log)),scale = 'row',main='heatmap')

write.csv(selY,file="select_gene_raw.csv",row.names = T)
write.csv(cpm_log2_sel,file="select_gene_log2cpm.csv",row.names = T)
write.csv(selected_list[[1]],file = "DEG_namelist.csv",row.names = F)



