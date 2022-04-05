setwd("D:/Research/Data")

####Read DEG list and correlation result data frame
DEG_list <-  read.csv("D:/Research/Data/DEG_namelist.csv", header = TRUE)
O_correlation <- read.csv("D:/Research/Data/Redo/Correlation_SCC_origin_rna.csv", header = TRUE,row.names = 1)

####seperate the DEG and non-DEG
DEG_data <- O_correlation[which(row.names(O_correlation)%in%DEG_list$x),]
non_DEG_data <- O_correlation[-which(row.names(O_correlation)%in%DEG_list$x),]

####Select the specific genus list(by DEG SCC>=0.3)
specific_genus <- c()
count <- 1
for (i in 1:length(names(DEG_data))) {
  if((max(DEG_data[names(DEG_data)[i]])>=0.3) | (min(DEG_data[names(DEG_data)[i]])<=(-0.3))){
    specific_genus[count] <- names(DEG_data)[i]
    count <- count+1
  }
}

####Multiple Fisher exact test
#count1 = abs(SCC)>=0.3 & deg
#count2 = abs(SCC)>=0.3 & non_deg
#count3 = abs(SCC)<0.3 & deg
#count4 = abs(SCC)<0.3 & non_deg
p_value_list<- c()
odd_ratio <- c()
for (i in 1:length(specific_genus)) {
  count1 <- 0
  count3 <- 0
  for (j in 1:length(rownames(DEG_data[specific_genus[i]]))) {
    if(abs(DEG_data[j,specific_genus[i]]) >= 0.3){
      count1 <- count1+1
    }
    else{
      count3 <- count3+1
    }
  }
  count2 <- 0
  count4 <- 0
  for (j in 1:length(rownames(non_DEG_data[specific_genus[i]]))) {
    if(abs(non_DEG_data[j,specific_genus[i]]) >= 0.3){
      count2 <- count2+1
    }
    else{
      count4 <- count4+1
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

combind <- cbind(p,specific_genus,odd_ratio)
combind <- data.frame(combind)
combind$p <- as.numeric(combind$p)
combind <- combind[which((combind$p < 0.05)&(combind$odd_ratio >1)),]
combind <- combind[order(combind$specific_genus),]
names(combind)<-c("q-value","genus","odds_ratio")

------------------
#compare with WGS results
wgs_result <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/gene_selected_genus_cutoff_2.csv",header = T)

library(ggvenn)
x<- list(WGS_results = wgs_result$Genus,
         RNA_results = combind$genus)

ggvenn(
  x, 
  fill_color = c("#FFCC99", "#99CCFF"),
  stroke_color = "#C0C0C0",
  stroke_size = 1, set_name_size = 8,
  label_sep = ",",
  text_size = 5
)
-------------------------------------
  
# Selected the odds ratio > 1 genus to rank DEG
selected_genus <- c(4,7,12,13,14,15,16,18,20,21,22,24,25,26,27,29,30,32,38,42,43,46,47)
selected_genus <- combind$genus[selected_genus]
select_DEG_table <- DEG_data[,which(names(DEG_data)%in%selected_genus)]
name_r <- rownames(select_DEG_table)
setwd("D:/Research/Data/GSEA_prerank")
for (i in 1:23) {
  table <-name_r
  table <- data.frame(table)
  table$scc <- select_DEG_table[,i]
  name <- names(select_DEG_table)[i]
  
  names(table)[1] <- name
  table <- table[order(-table$scc),]
  name <- strsplit(name,split = "\\.")
  name <- name[[1]][length(name[[1]])]
  filename <- paste(name,".rnk",sep = "")
  
  write.table(table,file=filename,quote=F,sep="\t",row.names=F,col.names = F)
}

