setwd("D:/Research/Data")

####Read DEG list and correlation result data frame
DEG_list <-  read.csv("D:/Research/Data/DEG_namelist.csv", header = TRUE)
O_correlation <- read.csv("D:/Research/Data/Redo/Correlation_SCC_origin.csv", header = TRUE,row.names = 1)

####seperate the DEG and non-DEG
DEG_data <- O_correlation[which(row.names(O_correlation)%in%DEG_list$x),]
non_DEG_data <- O_correlation[-which(row.names(O_correlation)%in%DEG_list$x),]

#####Calculate the count of 4 value
#count1 = abs(SCC)>=0.3 & deg
#count2 = abs(SCC)>=0.3 & non_deg
#count3 = abs(SCC)<0.3 & deg
#count4 = abs(SCC)<0.3 & non_deg
count1 <- 0
count3 <- 0
for (i in 1:length(names(DEG_data))) {
  for (j in 1:length(row.names(DEG_data))) {
    if(DEG_data[[i]][j] <= (-0.3)){
      count1 <- count1+1
    }
    else{
      count3 <- count3+1
    }
  }
}

count2 <- 0
count4 <- 0
for (i in 1:length(names(non_DEG_data))) {
  for (j in 1:length(row.names(non_DEG_data))) {
    if(non_DEG_data[[i]][j] <= (-0.3)){
      count2 <- count2+1
    }
    else{
      count4 <- count4+1
    }
  }
}

#combine 4 count together
counter <- rbind(c(count1,count2),
                 c(count3,count4))


result = fisher.test(counter)
