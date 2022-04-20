setwd("D:/Research/Data/Redo")
gene_select <- read.csv("Selected_genus_in_category/gene_selected_genus_cutoff_2.csv",header = T)
SCC_table <- read.csv("Correlation_SCC_origin_rna.csv",header = T,row.names = 1)
DEG_list<- read.csv("D:/Research/Data/DEG_namelist.csv",header = T)
SCC_select <- SCC_table[,which(names(SCC_table)%in%gene_select$Genus)]

####Multiple Fisher exact test
#count1 = abs(SCC)>=0.3 & deg
#count2 = abs(SCC)>=0.3 & non_deg
#count3 = abs(SCC)<0.3 & deg
#count4 = abs(SCC)<0.3 & non_deg
p_value_list<- c()
odd_ratio <- c()
for (i in 1:length(gene_select$Genus)) {
  count1 <- 0
  count2 <- 0
  count3 <- 0
  count4 <- 0
  for (j in 1:length(rownames(SCC_select))) {
    if(rownames(SCC_select)[j]%in%DEG_list$x){
      if(abs(SCC_select[j,gene_select$Genus[i]])>=0.3){
        count1<-count1+1
      }
      else{
        count3 <- count3+1
      }
    }
    else{
      if(abs(SCC_select[j,gene_select$Genus[i]])>=0.3){
        count2 <- count2+1
      }
      else{
        count4 <- count4+1
      }
    }
  }
  counter <- rbind(c(count1,count2),
                   c(count3,count4))
  
  
  result = fisher.test(counter)
  p_value_list[i] <- result$p.value
  odd_ratio[i] <- result$estimate
}

p = p_value_list
p = p.adjust(p,method="fdr",length(p)) 

combind <- cbind(gene_select$Genus,p,odd_ratio)
combind <- data.frame(combind)
combind$p <- as.numeric(combind$p)
names(combind)<-c("genus","q-value","odds_ratio")
