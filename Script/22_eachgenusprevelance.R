setwd("D:/Research/Data")
raw_abundance <-  read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)

#Remove the contaminant column
raw_abundance <- raw_abundance[,1:1284]
#Calculate the mean of all genus abundance
all_genus_mean <- apply(raw_abundance,2,mean)
#Plot distribution
hist(all_genus_mean)

#Whole matrix value
library(tibble)
library(tidyr)
all_genus_value <- rownames_to_column(raw_abundance, var = "Xval")
all_genus_value <- all_genus_value %>% pivot_longer(cols = 2:1285, names_to = "Source", values_to = "Value")

bin <- seq(-10,25)
hist(all_genus_value$Value,main = "Relative abundace values",xlab = "Values (Relative abundance)",ylab = "counts",breaks = bin)
abline(v = 2, col = "red", lwd = 3)

#Count the prevalence of each genus
Prevalence_genus <- function(cutoff){
  #The cutoff is the setting how much value is non-appear in human body
  result <- data.frame(names(raw_abundance))
  result$prevalence <- 0
  for (i in 1:length(names(raw_abundance))) {
    count_have_genus <- sum(raw_abundance[i]>cutoff) #count the number of sample that value > cutoff
    prevalence_temp <- count_have_genus/length(rownames(raw_abundance)) #calculate the prevalence
    result$prevalence[which(result$names.raw_abundance.==names(raw_abundance)[i])] <- prevalence_temp #fill the prevalence of the genus into reult
  }
  return(result)
}

zero <- Prevalence_genus(0)
one <- Prevalence_genus(1)
two <- Prevalence_genus(2)
test <-Prevalence_genus(3.5)
hist(test$prevalence,main = "Histogram of prevalence distribution",xlab = "Prevalence (%)",ylab = "Counts")
#plot the prevalence of genus
library(ggplot2)
test <- test[order(test$prevalence),]
test$num <- c(1:1284)
test$prevalence_exceed_20 <- ifelse(test$prevalence>0.2,T,F)
ggplot(data = test,aes(x=num,y=prevalence)) +
  geom_line(size = 1)+
  geom_hline(yintercept =0.2,color = "red")+
  labs(title = "Prevalence of each genus",subtitle = "Cutoff=3.5",x="Genus",y="Prevalence")

sum(test$prevalence_exceed_20 == T)
write.csv(zero,file="D:/Research/Data/Redo/Cutoff_0_genus.csv",row.names = F)
write.csv(one,file="D:/Research/Data/Redo/Cutoff_1_genus.csv",row.names = F)
write.csv(two,file="D:/Research/Data/Redo/Cutoff_2_genus.csv",row.names = F)
