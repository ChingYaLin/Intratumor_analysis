setwd("D:/Research/Data/Redo")

miRNA_bacter <- read.csv("Bacteria_miRNA.csv",row.names = 1)
my_bacter <- read.csv("Bacteria_mydata.csv",row.names = 1)

test <- my_bacter[which(my_bacter$phylum!=""),]
miRNA_raw <- read.csv("D:/Research/Data/Raw/miRNA microbiome/LUAD_relative_abundance_miRNA.txt",sep = "\t")
nor_data <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)
raw_data <- read.csv("D:/Research/Data/Raw/Kraken-TCGA-Raw-Data-17625-Samples.csv", header = TRUE,row.names = 1)

samples <- rownames(nor_data)
genus <- names(nor_data)

raw_138 <- raw_data[which(rownames(raw_data)%in%samples),which(names(raw_data)%in%genus)]

#remove contaminant (last 3 column)
nor_data <- nor_data[,1:1284]
cut_genus <- nor_data



#Calculate the correlation between miRNA and WGS genus
s_genus <- read.table("D:/Research/Data/Raw/miRNA microbiome/LUAD_relative_abundance_miRNA.txt",header = TRUE, sep="\t")
temp <- colnames(s_genus)
temp <- substr(temp,1,12)
temp <- gsub("[.]","-",temp)
names(s_genus) <- temp

tgca <- read.csv("D:/Research/TCGA/TCGA.csv", header = TRUE,row.names = 1)
s_genus <-s_genus[,which(names(s_genus)%in%tgca$barcode)]

for (i in 1:length(names(s_genus))) {
  names(s_genus)[i] <- rownames(tgca)[which(tgca$barcode==names(s_genus)[i])]
}
s_genus<-s_genus[,1:134]

raw_134 <- raw_138[which(rownames(raw_138)%in%names(s_genus)),]
#rename the genus name by genus level
colname_genus <- colnames(cut_genus)
for (j in 1:length(colname_genus)) {
  if (colname_genus[j]!='X'){
    temp <- strsplit(colname_genus[j],split = "\\.")
    names(cut_genus)[names(cut_genus)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
  }
}
temp <- gsub("g_","",names(cut_genus))
temp <- gsub("[_]","",temp)
names(cut_genus) <- temp

colname_genus <- colnames(raw_134)
for (j in 1:length(colname_genus)) {
  if (colname_genus[j]!='X'){
    temp <- strsplit(colname_genus[j],split = "\\.")
    names(raw_134)[names(raw_134)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
  }
}
temp <- gsub("g_","",names(raw_134))
temp <- gsub("[_]","",temp)
names(raw_134) <- temp

s_genus <-s_genus[which(rownames(s_genus)%in%colnames(cut_genus)),]
cut_genus <- cut_genus[,which(names(cut_genus)%in%rownames(s_genus))]
cut_genus <- cut_genus[which(rownames(cut_genus)%in%names(s_genus)),]
raw_134 <- raw_134[,which(names(raw_134)%in%rownames(cut_genus))]

####spearman correlation coefficient
#create a blank dataframe
s_genus <- data.frame(t(s_genus))
s_genus<- s_genus[ order(row.names(s_genus)), ]
cut_genus <- cut_genus[ order(row.names(cut_genus)), ]

s_genus<- s_genus[ , order(names(s_genus))]
cut_genus <- cut_genus[ , order(names(cut_genus))]

s_genus <- data.frame(t(s_genus))
cut_genus <- data.frame(t(cut_genus))

result <- names(cut_genus)
result <- data.frame(result)
result$pvalue <- 0
result$SCC <- 0

#calculate the correlation and fill the matrix
for (i in 1:length(names(s_genus))) {
  m<-names(s_genus)[i]
  g<-names(cut_genus)[i]
  test<-cor.test(s_genus[[m]], cut_genus[[g]], method=c("spearman"),exact=F)
  rho <- test$estimate[["rho"]]
  pvalue <- test$p.value
  result$SCC[i] <- rho
  result$pvalue[i] <- pvalue
}
hist(result$SCC)
test <- na.omit(result)


####spearman correlation coefficient
raw_134 <- data.frame(t(raw_134))
result <- names(raw_134)
result <- data.frame(result)
result$pvalue <- 0
result$SCC <- 0

#calculate the correlation and fill the matrix
for (i in 1:length(names(s_genus))) {
  m<-names(s_genus)[i]
  g<-names(raw_134)[i]
  test<-cor.test(s_genus[[m]], raw_134[[g]], method=c("spearman"),exact=F)
  rho <- test$estimate[["rho"]]
  pvalue <- test$p.value
  result$SCC[i] <- rho
  result$pvalue[i] <- pvalue
}
hist(result$SCC)
test <- na.omit(result)


s_genus$phylum <- miRNA_bacter$phylum[which(miRNA_bacter$query%in%rownames(s_genus))]
raw_134$phylum <- my_bacter$phylum[which(my_bacter$query%in%rownames(raw_134))]
library(data.table)

DT1 <- data.table(s_genus)  # DF is your original data
library(plyr)
temp <- ddply(s_genus, "phylum", numcolwise(sum))
temp2 <- ddply(raw_134,"phylum",numcolwise(sum))

temp2 <- temp2[1:35,]

rownames(temp) <- temp$phylum
rownames(temp2) <- temp2$phylum
temp <- temp[,-1]
temp2 <- temp2[,-1]
temp2 <- temp2[,order(names(temp2))]
temp_origin <- data.frame(t(temp2))
temp_origin <- temp_origin/rowSums(temp_origin)
temp_origin <- data.frame(t(temp_origin))

temp_miRNA <- data.frame(t(temp))
temp_miRNA <- temp_miRNA/rowSums(temp_miRNA)
temp_miRNA <- data.frame(t(temp_miRNA))
temp_miRNA <- temp_miRNA[,-which(names(temp_miRNA)=="s11217")]
temp_origin <- temp_origin[,-which(names(temp_origin)=="s11217")]

##phylum correlation(by sample)
result <- names(temp_miRNA)
result <- data.frame(result)
result$pvalue <- 0
result$SCC <- 0

#calculate the correlation and fill the matrix
for (i in 1:length(names(temp_miRNA))) {
  m<-names(temp_miRNA)[i]
  g<-names(temp_origin)[i]
  test<-cor.test(temp_miRNA[[m]], temp_origin[[g]], method=c("spearman"),exact=F)
  rho <- test$estimate[["rho"]]
  pvalue <- test$p.value
  result$SCC[i] <- rho
  result$pvalue[i] <- pvalue
}
hist(result$SCC,main = "Compare correlation between miRNA and WGS sample",
     xlab = "SCC",
     ylab = "frequency")
test <- na.omit(result)

##phylum correlation(by phylum)
temp_miRNA_g <- data.frame(t(temp_miRNA))
temp_origin_g <- data.frame(t(temp_origin))
result <- names(temp_miRNA_g)
result <- data.frame(result)
result$pvalue <- 0
result$SCC <- 0

#calculate the correlation and fill the matrix
for (i in 1:length(names(temp_miRNA_g))) {
  m<-names(temp_miRNA_g)[i]
  g<-names(temp_origin_g)[i]
  test<-cor.test(temp_miRNA_g[[m]], temp_origin_g[[g]], method=c("spearman"),exact=F)
  rho <- test$estimate[["rho"]]
  pvalue <- test$p.value
  result$SCC[i] <- rho
  result$pvalue[i] <- pvalue
}
hist(result$SCC,main = "Compare correlation between miRNA and WGS sample",
     xlab = "SCC",
     ylab = "frequency")
test <- na.omit(result)
