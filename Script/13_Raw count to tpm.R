library (EDASeq)
library("edgeR")
setwd("D:/Research/Data")
#load data
raw_count <-  read.csv("D:/Research/Data/Raw/LUAD_gene_expression_raw(138sample).csv", header = TRUE)
#remove normal sample
raw_count <- raw_count[,-c(140:173)]
raw_count <- raw_count[-c(60484:60488),]
#Rename the column
name <- names(raw_count)
for (i in 1:length(name)) {
  temp_len <- nchar(name[i])
  temp <- substr(name[i],temp_len-5,temp_len)
  names(raw_count)[names(raw_count) == name[i]] <- temp
}
#Rename the raw name
rownames(raw_count) <- raw_count$X
raw_count <- raw_count[,-which(names(raw_count)%in%("X"))]

cpm_log <- cpm(raw_count, log = TRUE)
mean_log2_cpm <- apply(cpm_log, 1, mean)
hist(mean_log2_cpm)
expr_cutoff <- 0
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(mean_log2_cpm > expr_cutoff)
data_clean <- raw_count[mean_log2_cpm > expr_cutoff, ]

ensembl_list <- rownames(data_clean)
ensembl_list<- ensembl_list[15001:15329]
result <- getGeneLengthAndGCContent(ensembl_list, "hsa")

#resultt <- read.csv(file="D:/Research/Data/Redo/Gene_length.csv",row.names = 1)
temp <- getGeneLengthAndGCContent(ensembl_list, "hsa")
resultt <- rbind(resultt, temp)
resultt <- data.frame(resultt)
#remove na rows
result <- na.omit(resultt) 
#remove the rows that delete by na
test<-data_clean
test <- test[which(rownames(test)%in%rownames(result)),]
#Convert ensembl ID to HGNC gene symbol
library(biomaRt)
ensembl_list <- rownames(test)
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_coords=getBM(attributes=c("hgnc_symbol","ensembl_gene_id"), filters="ensembl_gene_id", values=ensembl_list, mart=human)
#Rename the rownames
test$X <- rownames(test)
test$HGNC<-0
for (i in 1:length(test$HGNC)) {
  test$HGNC[i]<-gene_coords[which(gene_coords$ensembl_gene_id==test$X[i]),"hgnc_symbol"]
}
#Remove the rows that not have HGNC symbol
test <- test[-which(test$HGNC%in%c("")),]
result <- result[which(rownames(result)%in%test$X),]

test <- test[-which(rownames(test)=="ENSG00000272655"),]
result <- result[-which(rownames(result)=="ENSG00000272655"),]
rownames(test) <-test$HGNC
test <- test[,-which(names(test)%in%c("HGNC","X"))]
#convert to tpm
for (i in 1:length(names(test))) {
  countt<-test[names(test)[i]]
  len <- result$length
  test[names(test)[i]] <- countt/len
  test[names(test)[i]] <- test[names(test)[i]]/sum(test[names(test)[i]],na.rm = T)*1e6
}
write.csv(test,file="D:/Research/Data/Redo/138sample_tpm_HGNC1.csv",row.names = T)
write.csv(resultt,file="D:/Research/Data/Redo/Gene_length.csv",row.names = T)


