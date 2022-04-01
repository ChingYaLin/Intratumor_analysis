setwd("D:/Research/Data/Raw")
df_wgs <- read.csv("LUAD_WGS_abundance.csv")
df_rna <- read.csv("LUAD_RNA_seq_abundance.csv")
Colname <- colnames(df_wgs)

#rename the column's name to genus level
for (j in 1:length(Colname)) {
  if (Colname[j]!='X'){
    temp <- strsplit(Colname[j],split = "\\.")
    names(df_rna)[names(df_rna)==Colname[j]] <- temp[[1]][length(temp[[1]])]
    names(df_wgs)[names(df_wgs)==Colname[j]] <- temp[[1]][length(temp[[1]])]
  }
}


#remove the column of sample name and contaminant
drop <- c("X","contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters")
df_clear_rna = df_rna[,!(names(df_rna) %in% drop)]
df_clear_wgs = df_wgs[,!(names(df_wgs) %in% drop)]

# calculate the prevalence of each bacteria
Colname <- colnames(df_clear_rna)
list_rna <- c()
for (i in 1:length(Colname)) {
  count <- 0
  for (j in 1:519) {
    if(df_clear_rna[j,Colname[i]]>0){
      count <- count+1
    }
  }
  prevalence <- count/519
  list_rna[i]<-prevalence
}
list_wgs <- c()
for (i in 1:length(Colname)) {
  count <- 0
  for (j in 1:162) {
    if(df_clear_wgs[j,Colname[i]]>0){
      count <- count+1
    }
  }
  prevalence <- count/162
  list_wgs[i]<-prevalence
}

#draw the prevalence curve
library(ggplot2)
num <- c(1:length(Colname))
data <- rbind(list_rna, list_wgs)
data <- data.frame(data)
for (i in 1:length(colnames(data))) {
  names(data)[names(data) == colnames(data)[i]] <- Colname[i]
}
data <- t(data)
data <- data.frame(data)
newdata <- data[order(-list_rna,-list_wgs),]
newdata <- data
newdata <- t(newdata)
newdata<-rbind(newdata,num)
newdata <- data.frame(t(newdata))


ggplot(newdata) +
  geom_line(aes(num, list_rna),col = "red")+
  geom_line(aes(num, list_wgs), col = "blue")+
  labs(title = "Prevalence", subtitle = "red line: RNA-seq / blue line: WGS",x="Microbiome", y="prevalence")


##############################################################################################################
#same sample
setwd("D:/Research")
df_wgs <- read.csv("LUAD_WGS_abundance_only138.csv")
df_rna <- read.csv("LUAD_RNA_seq_abundance_only138.csv")
Colname <- colnames(df_wgs)

#rename the column's name to genus level
for (j in 1:length(Colname)) {
  if (Colname[j]!='X'){
    temp <- strsplit(Colname[j],split = "\\.")
    names(df_rna)[names(df_rna)==Colname[j]] <- temp[[1]][length(temp[[1]])]
    names(df_wgs)[names(df_wgs)==Colname[j]] <- temp[[1]][length(temp[[1]])]
  }
}


#remove the column of sample name and contaminant
drop <- c("X","contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters")
df_clear_rna = df_rna[,!(names(df_rna) %in% drop)]
df_clear_wgs = df_wgs[,!(names(df_wgs) %in% drop)]

# calculate the prevalence of each bacteria
Colname <- colnames(df_clear_rna)
list_rna <- c()
for (i in 1:length(Colname)) {
  count <- 0
  for (j in 1:138) {
    if(df_clear_rna[j,Colname[i]]>0){
      count <- count+1
    }
  }
  prevalence <- count/138
  list_rna[i]<-prevalence
}
list_wgs <- c()
for (i in 1:length(Colname)) {
  count <- 0
  for (j in 1:138) {
    if(df_clear_wgs[j,Colname[i]]>0){
      count <- count+1
    }
  }
  prevalence <- count/138
  list_wgs[i]<-prevalence
}

#draw the prevalence curve
library(ggplot2)
num <- c(1:length(Colname))
data <- rbind(list_rna, list_wgs)
data <- data.frame(data)
for (i in 1:length(colnames(data))) {
  names(data)[names(data) == colnames(data)[i]] <- Colname[i]
}
data <- t(data)
data <- data.frame(data)
newdata <- data[order(-list_rna,-list_wgs),]
newdata <- t(newdata)
newdata<-rbind(newdata,num)
newdata <- data.frame(t(newdata))
newdata$compare <- newdata$list_rna-newdata$list_wgs


ggplot(newdata) +
  geom_histogram(aes(compare),fill = "blue",binwidth = 0.01)+
  labs(title = "Prevalence compare", subtitle = "RNA-seq prevalence - WGS prevalence",x="value", y="count")

ggplot(newdata,aes(x=num,y=compare))+
  geom_line(color="blue")+
  geom_hline(yintercept = 0,color = "red")+
  labs(title = "Prevalence compare", subtitle = "RNA-seq prevalence - WGS prevalence",x="microbiome", y="comparison")

#count the number of > and <
rna_bigger <- 0
wgs_bigger <- 0
equal <- 0
for (i in 1:length(newdata$compare)) {
  if(newdata$compare[i]>0.1){
    rna_bigger <- rna_bigger+1
  }
  else if(newdata$compare[i]<(-0.1)){
    wgs_bigger <- wgs_bigger+1
  }
  else{
    equal <- equal+1
  }
}
cat(rna_bigger,wgs_bigger,equal)

##############################################################################################################

