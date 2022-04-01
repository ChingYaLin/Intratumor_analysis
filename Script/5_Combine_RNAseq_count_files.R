setwd("D:/Research/TCGA")
library("dplyr")
manifast <- read.table("gdc_manifest_20210628_120505.txt",header = T)
allsample <- read.csv("D:/Research/TCGA/TCGA.csv",header = T)

tumor_dir <- allsample$tumor
normal_dir<- allsample$normal
normal_dir <- normal_dir[!normal_dir %in% c("")]
fil <- allsample %>% filter(normal %in% normal_dir)
tumor_c <- fil$tumor

#merge the tumor sample count together and rename columns
count<-2
test <- paste("D:/Research/TCGA",manifast$id[1],manifast$filename[1],sep = '/')
sample_name <- allsample$X[which(allsample$tumor==manifast$id[1])]
com <- read.table(test)
names(com)[names(com) == "V2"] <- paste("TC_1_",sample_name,sep = "")
for (i in 2:length(manifast$id)) {
  if(manifast$id[i] %in% tumor_dir){
    test <- paste("D:/Research/TCGA",manifast$id[i],manifast$filename[i],sep = '/')
    temp <- temp <- read.table(test)
    temp$V1<-substr(temp$V1,1,15)
    rownames(temp) <- temp$V1
    if(manifast$id[i] %in% tumor_c){
      sample_name <- allsample$X[which(allsample$tumor==manifast$id[i])]
      name <- paste("TC_",as.character(count),"_",sample_name,sep = "")
      names(temp)[names(temp) == "V2"] <- name
      count <- count+1
    }
    else{
      sample_name <- allsample$X[which(allsample$tumor==manifast$id[i])]
      name <- paste("TS_",as.character(count),"_",sample_name,sep = "")
      names(temp)[names(temp) == "V2"] <- name
      count <- count+1
    }
    com <- cbind(com,temp)
    com<-com[,!(names(com) %in% c("V1"))]
  }
}

#merge the normal sample file into com and rename by N_1,N_2...
count<-1
for (i in 1:length(manifast$id)) {
  if(manifast$id[i] %in% normal_dir){
    test <- paste("D:/Research/TCGA",manifast$id[i],manifast$filename[i],sep = '/')
    temp <- temp <- read.table(test)
    temp$V1<-substr(temp$V1,1,15)
    rownames(temp) <- temp$V1
    name <- paste("NO_",as.character(count),sep = "")
    names(temp)[names(temp) == "V2"] <- name
    count <- count+1
    com <- cbind(com,temp)
    com<-com[,!(names(com) %in% c("V1"))]
  }
}

write.csv(com,file = 'D:/Research/Data/LUAD_DEG_raw.csv',row.names = T)
