setwd("D:/Research/Data/Redo")
nor_data <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)
raw_data <- read.csv("D:/Research/Data/Raw/Kraken-TCGA-Raw-Data-17625-Samples.csv", header = TRUE,row.names = 1)

samples <- rownames(nor_data)
genus <- names(nor_data)

raw_138 <- raw_data[which(rownames(raw_data)%in%samples),which(names(raw_data)%in%genus)]

#remove contaminant (last 3 column)
nor_data <- nor_data[,1:1284]


min(nor_data)
mean_normalized_data <- apply(nor_data, 2, mean)
sum(mean_normalized_data > 0)
cut_genus <- nor_data[, mean_normalized_data > 0]



#Calculate the correlation between miRNA and WGS genus
s_genus <- read.table("D:/Research/Data/Raw/LUAD_relative_abundance_miRNA.txt",header = TRUE, sep="\t")
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

s_genus <-s_genus[which(rownames(s_genus)%in%colnames(cut_genus)),]
cut_genus <- cut_genus[,which(names(cut_genus)%in%rownames(s_genus))]
cut_genus <- cut_genus[which(rownames(cut_genus)%in%names(s_genus)),]

####spearman correlation coefficient
#create a blank dataframe
s_genus <- data.frame(t(s_genus))
s_genus<- s_genus[ order(row.names(s_genus)), ]
cut_genus <- cut_genus[ order(row.names(cut_genus)), ]

s_genus<- s_genus[ , order(names(s_genus))]
cut_genus <- cut_genus[ , order(names(cut_genus))]

s_genus <- data.frame(t(s_genus))

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


