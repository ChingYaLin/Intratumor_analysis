#### Set pathway ####
dir = "D:/Research/TCGA-LUAD all case"
setwd(dir)

#### Load data####
all_luad <- read.csv("gdc_sample_sheet.2022-02-09.tsv", header = T,sep = "\t")
my_luad <- read.csv("D:/Research/TCGA/TCGA.csv",header = T)

#### Extract Normal sample barcode and check####
all_normal <- all_luad[which(all_luad$Sample.Type == "Solid Tissue Normal"),]
all_tumor <- all_luad[which(all_luad$Sample.Type != "Solid Tissue Normal"),]
match_sample <- all_normal$Case.ID
remain <- match(match_sample, all_tumor$Case.ID)
remain_tumor <- all_tumor[remain,]
my_nor <- my_luad[which(my_luad$normal!=""),]
my_nor_code <- my_nor$barcode
all_nor_code <-all_normal$Case.ID
my_nor_code<- data.frame(my_nor_code)
my_nor_code$exsist <- ifelse(my_nor_code$my_nor_code %in% remain_tumor$Case.ID,TRUE,FALSE)

####Combine RNAseq count files####
#merge tumor sample files and rename by T_1, T_2, ...
read_dir <- paste("D:/Research/TCGA-LUAD all case",all_tumor$File.ID[1],all_tumor$File.Name[1],sep = '/')
count<-2
com <- read.table(read_dir)
names(com)[names(com) == "V2"] <- "T_1"
com$V1 <- substr(com$V1,1,15)
rownames(com) <- com$V1
for (i in 2:length(all_tumor$File.ID)) {
  read_dir <- paste("D:/Research/TCGA-LUAD all case",all_tumor$File.ID[i],all_tumor$File.Name[i],sep = '/')
  temp <- read.table(read_dir)
  rownames(temp)<-substr(temp$V1,1,15)
  name <- paste("T_",count,sep = "")
  names(temp)[names(temp) == "V2"] <- name
  count <- count+1
  com <- cbind(com,temp)
  com<-com[,!(names(com) %in% c("V1"))]
}
#merge the normal sample file into com and rename by N_1,N_2...
count<-1
for (i in 1:length(all_normal$File.ID)) {
  read_dir <- paste("D:/Research/TCGA-LUAD all case",all_normal$File.ID[i],all_normal$File.Name[i],sep = '/')
  temp <- read.table(read_dir)
  rownames(temp) <- substr(temp$V1,1,15)
  name <- paste("N_",count,sep = "")
  names(temp)[names(temp) == "V2"] <- name
  count <- count+1
  com <- cbind(com,temp)
}
com<-com[,!(names(com) %in% c("V1"))]
write.csv(com,file = 'D:/Research/TCGA-LUAD all case/LUAD_DEG_raw.csv',row.names = T)
