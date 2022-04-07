setwd("D:/Research/Data")

#Read metadata that have live and dead informations
meta <-  read.csv("D:/Research/Data/Raw/Metadata_LUAD_WGS_only138.csv", header = TRUE)
#Seperate the alive and dead group
alive <- meta[which(meta$vital_status=="Alive"),]
dead <- meta[which(meta$vital_status=="Dead"),]

alive <- alive$X
dead <- dead$X
#Load the genus abundance matrix
genus <- read.csv("D:/Research/Data/Redo/voom_svm_clean_genus.csv", header = TRUE,row.names = 1)

alive_g <- genus[which(rownames(genus)%in%alive),]
dead_g <- genus[which(rownames(genus)%in%dead),]

ttest_list <- c()
for (i in 1:length(names(genus))) {
 ttest<- t.test(alive_g[names(genus)[i]],dead_g[names(genus)[i]]) 
 ttest_list[i] <- ttest$p.value
}

ttest_list <- data.frame(ttest_list)
rownames(ttest_list) <- names(genus)

count <- 0
for (i in 1:length(ttest_list$ttest_list)) {
  if(ttest_list$ttest_list[i]<0.05){
    count <- count+1
  }
}

ttest_list$genus <- rownames(ttest_list)
significant_genus <- ttest_list[which(ttest_list$ttest_list<0.05),]

p <- ttest_list$ttest_list
p<-p.adjust(p,method="fdr",length(p))
p <- data.frame(p)
p$genus<-ttest_list$genus

test<-p[which(p$p<0.05),]
