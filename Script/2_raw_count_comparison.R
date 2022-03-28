library(tidyverse)
setwd("D:/Research/Data")
df_raw <- read.csv("Kraken-TCGA-Raw-Data-17625-Samples.csv")
filter_rna <- read.csv("Metadata_LUAD_RNA_only138.csv")
filter_wgs <- read.csv("Metadata_LUAD_WGS_only138.csv")
df_kraken <- read.csv("Kraken-TCGA-Voom-SNM-Plate-Center-Filtering-Data.csv")
df_meta <- read.csv("Metadata-TCGA-Kraken-17625-Samples.csv")

#select 1285 columns form 1994 columns
drop <- c("contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters")
df_kraken = df_kraken[,!(names(df_kraken) %in% drop)]
colname <- colnames(df_kraken)
df <- df_raw %>% select(colname)
df$UUID <- df_meta$sample_uuid

sample_list_RNA<- filter_rna$X
data_filter_rna <- df[df$X %in% c(sample_list_RNA),]
sample_list_WGS<- filter_wgs$X
data_filter_wgs <- df[df$X %in% c(sample_list_WGS),]



#rename the bacteria name into genus
colname_rna <- colnames(data_filter_rna)
colname_wgs <- colnames(data_filter_wgs)
for (j in 1:length(colname_rna)) {
  if (colname_rna[j]!='X'){
    temp <- strsplit(colname_rna[j],split = "\\.")
    names(data_filter_rna)[names(data_filter_rna)==colname_rna[j]] <- temp[[1]][length(temp[[1]])]
  }
}
for (j in 1:length(colname_wgs)) {
  if (colname_wgs[j]!='X'){
    temp <- strsplit(colname_wgs[j],split = "\\.")
    names(data_filter_wgs)[names(data_filter_wgs)==colname_wgs[j]] <- temp[[1]][length(temp[[1]])]
  }
}

#remove the X column
data_filter_rna <- data_filter_rna[,!(names(data_filter_rna) %in% c("X"))]
data_filter_wgs <- data_filter_wgs[,!(names(data_filter_wgs) %in% c("X"))]

#sum the row by sample
rna_sum <- rowSums(data_filter_rna[,1:1284])
rna_sum <- data.frame(rna_sum)
rna_sum$X <- data_filter_rna$UUID

wgs_sum <- rowSums(data_filter_wgs[,1:1284])
wgs_sum <- data.frame(wgs_sum)
wgs_sum$X <- data_filter_wgs$UUID

rna_sum<-rna_sum[str_order(rna_sum$X),]
wgs_sum<-wgs_sum[str_order(wgs_sum$X),]

rownames(rna_sum) <- rna_sum$X
rownames(wgs_sum) <- wgs_sum$X
rna_sum <- rna_sum[,!(names(rna_sum) %in% c("X"))]
wgs_sum <- wgs_sum[,!(names(wgs_sum) %in% c("X"))]
test <- rbind(wgs_sum,rna_sum)
test <- t(test)
test<-data.frame(test)

ggplot(test,aes(x=wgs_sum,y=rna_sum))+
  geom_point(colour = "blue", size = 2,shape = 21)+
  geom_hline(yintercept = 500000,color = "red")+
  geom_vline(xintercept = 500000,color="red")+
  labs(title = "Raw counts compare",x="WGS", y="RNA-seq")


#count the number of > and <
ttest <- test$wgs_sum/test$rna_sum
rna_bigger <- 0
wgs_bigger <- 0
equal <- 0
for (i in 1:length(ttest)) {
  if(ttest[i]<=0.5){
    rna_bigger <- rna_bigger+1
  }
  else if (ttest[i]>=2){
    wgs_bigger <- wgs_bigger+1
  }
  else{
    equal <- equal+1
  }
}
cat(rna_bigger,wgs_bigger,equal)

wgs_mean <- mean(test$wgs_sum)
rna_mean <- mean(test$rna_sum)
cat("WGS mean:",wgs_mean,"/ RNA-seq mean:",rna_mean)