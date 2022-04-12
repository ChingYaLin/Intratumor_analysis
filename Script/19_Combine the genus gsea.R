setwd("D:/Research/Data/GSEA_result")
genus_list <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/gene_selected_genus_cutoff_2.csv",header = T)
genus_list <- genus_list$Genus
#rename the genus name
for (j in 1:length(genus_list)) {
  temp <- strsplit(genus_list[j],split = "\\.")
  genus_list[j] <- temp[[1]][length(temp[[1]])]
}


neg<-read.table("gsea_report_for_na_neg_g__Actinotalea.tsv",sep = "\t")
pos<-read.table("gsea_report_for_na_pos_g__Actinotalea.tsv",sep = "\t")
#The last column is NA so need to remove
neg <- neg[,1:11]
pos <- pos[,1:11]
#The first row is column name so rename columns.
names(neg) <- neg[1,]
names(pos) <- pos[1,]
#Remove the first row that is name.
neg <- neg[-1,]
pos <- pos[-1,]
#Add the column that indicate which row is what kind of genus.
genus <- "g__Actinotalea"
neg$genus <- genus
pos$genus <- genus

final_combine_table_neg <- neg
final_combine_table_pos <- pos

####Start to deal with 23 genus's gsea results and combine
for (i in 2:length(genus_list)) {
  negpath <- paste("D:/Research/Data/GSEA_result/gsea_report_for_na_neg_",genus_list[i],".tsv",sep = "")
  pospath <- paste("D:/Research/Data/GSEA_result/gsea_report_for_na_pos_",genus_list[i],".tsv",sep = "")
  neg<-read.table(negpath,sep = "\t")
  pos<-read.table(pospath,sep = "\t")
  neg <- neg[,1:11]
  pos <- pos[,1:11]
  names(neg) <- neg[1,]
  names(pos) <- pos[1,]
  neg <- neg[-1,]
  pos <- pos[-1,]
  neg$genus <- genus_list[i]
  pos$genus <- genus_list[i]
  
  #Combine the dataframe
  final_combine_table_neg <- rbind(final_combine_table_neg,neg)
  final_combine_table_pos <- rbind(final_combine_table_pos,pos)
}



all_table <- rbind(final_combine_table_neg,final_combine_table_pos)
write.csv(final_combine_table_neg,file = "D:/Research/Data/Redo/GSEA_23_neg.csv")
write.csv(final_combine_table_pos,file = "D:/Research/Data/Redo/GSEA_23_pos.csv")
#create a blank dataframe
final_combine_table_pos <- read.csv("D:/Research/Data/Redo/GSEA_23_pos.csv", header = TRUE, row.names = 1)
final_combine_table_neg <- read.csv("D:/Research/Data/Redo/GSEA_23_neg.csv", header = TRUE, row.names = 1)
all_table <- rbind(final_combine_table_neg,final_combine_table_pos)
function_anno <- all_table$NAME
function_anno <- unique(function_anno)
genus <- final_combine_table_neg$genus
genus <- unique(genus)

heat_map <- data.frame(function_anno)
fill <- rep(0L, times = length(function_anno))
for (i in genus) {
  heat_map$x <- fill
  names(heat_map)[names(heat_map) == "x"] <- i
}
rownames(heat_map)<- heat_map$function_anno
heat_map = heat_map[,!(names(heat_map) %in% c("function_anno"))]

final_combine_table_neg$trans <- (-log10(as.numeric(final_combine_table_neg$FDR.q.val)))
final_combine_table_pos$trans <- (-log10(as.numeric(final_combine_table_pos$FDR.q.val)))

for (i in 1:length(final_combine_table_neg$trans)) {
  if(final_combine_table_neg$trans[i]==Inf){
    final_combine_table_neg$trans[i] <- 5
  }
}
for (i in 1:length(final_combine_table_pos$trans)) {
  if(final_combine_table_pos$trans[i]==Inf){
    final_combine_table_pos$trans[i] <- 5
  }
}
final_combine_table_neg$trans <- (-final_combine_table_neg$trans)
for (i in 1:length(final_combine_table_neg$NAME)) {
  heat_map[final_combine_table_neg$NAME[i],final_combine_table_neg$genus[i]] <-final_combine_table_neg$trans[i]
}
for (i in 1:length(final_combine_table_pos$NAME)) {
  heat_map[final_combine_table_pos$NAME[i],final_combine_table_pos$genus[i]] <-final_combine_table_pos$trans[i]
}

heatmap(t(t(heat_map)),scale="row")

library("pheatmap")
draw <- t(t(heat_map))
pheatmap(draw,cutree_rows = 4,show_colnames = F,fontsize_row = 3)
out <- pheatmap(draw,cutree_rows = 4,show_colnames = F,fontsize_row = 3,clustering_distance_cols="euclidean")
rownames(draw[out$tree_row[["order"]],])
colnames(draw[,out$tree_col[["order"]]])
sort(cutree(out$tree_row, k=4))
plot(out$tree_col)
abline(h=60, col="red", lty=2, lwd=2)
cluster = which(sort(cutree(out$tree_row, k=4))==3)
cluster_rm <- names(cluster)
cut_heat_map <- heat_map[-which(rownames(heat_map)%in%cluster_rm),]
pdf(file = "D:/Research/Plot/19.FunctionAnnotationAndGenus.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 30)
par(mar=c(10,10,3,3))
pheatmap(cut_heat_map,cutree_rows = 3,fontsize_row =3 ,show_colnames = T,
         main = "Genus associated functions")
dev.off()
cluster1 <- heat_map[which(rownames(heat_map)%in%cluster_rm),]
cluster2 <- heat_map[which(rownames(heat_map)%in%cluster_rm),]
cluster3 <- heat_map[which(rownames(heat_map)%in%cluster_rm),]
cluster4 <- heat_map[which(rownames(heat_map)%in%cluster_rm),]
pheatmap(cluster4,fontsize_row = 5)
