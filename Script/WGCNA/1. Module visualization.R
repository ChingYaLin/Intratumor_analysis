workingDir = "D:/Research/Script/WGCNA";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

####Load Data####
# Load the expression and trait data saved in the first part
lnames = load(file = "FemaleLiver-01-dataInput.RData");
# Load network data saved in the second part.
lnames = load(file = "FemaleLiver-02-networkConstruction-stepByStep.RData");
lnames = load(file = "Genus-01-dataInput.RData");
lnames = load(file = "Genus-02-networkConstruction-stepByStep.RData");


####Module correlation####
#Keep the same sample
remain <- match(rownames(MEs),rownames(MEs_g))
MEs_g = MEs_g[remain, ]
#remove the grey module
MEs <- MEs[,-which(names(MEs)%in%c("MEgrey"))]
MEs_g <- MEs_g[,-which(names(MEs_g)%in%c("MEgrey"))]
#Rename the genus module to add g at the module name
names(MEs_g) <- paste(names(MEs_g),"_g",sep="")
#create a blank dataframe
genusm_name <- colnames(MEs_g)
genem_name <- colnames(MEs)
correlation_module <- data.frame(genusm_name)
pvalue_module <- data.frame(genusm_name)
fill <- rep(0L, times = length(correlation_module$genusm_name))


for (i in genem_name) {
  correlation_module$x <- fill
  pvalue_module$x<- fill
  names(correlation_module)[names(correlation_module) == "x"] <- i
  names(pvalue_module)[names(pvalue_module) == "x"] <- i
}
rownames(correlation_module)<- correlation_module$genusm_name
correlation_module = correlation_module[,!(names(correlation_module) %in% c("genusm_name"))]
rownames(pvalue_module)<- pvalue_module$genusm_name
pvalue_module = pvalue_module[,!(names(pvalue_module) %in% c("genusm_name"))]
#calculate the correlation and fill the matrix
for (i in 1:length(genusm_name)) {
  for (j in 1:length(genem_name)) {
    m<-genusm_name[i]
    g<-genem_name[j]
    test<-cor.test(MEs[[g]], MEs_g[[m]], method=c("spearman"),exact=F)
    rho <- test$estimate[["rho"]]
    pvalue <- test$p.value
    correlation_module[which(rownames(correlation_module)==m),which(names(correlation_module)==g)] <- rho
    pvalue_module[which(rownames(pvalue_module)==m),which(names(pvalue_module)==g)] <- pvalue
  }
}

library(pheatmap)
pdf(file = "C:/Users/user/Desktop/My Plot1.pdf",   # The directory you want to save the file in
    width = 20, # The width of the plot in inches
    height = 10)
par(mar = c(6, 8.5, 3, 3))
pheatmap(correlation_module,main = "The correlation of genus Communities and gene modules")
dev.off()

####Module x module network####
library(tibble)
library(tidyr)
DF <- rownames_to_column(correlation_module, var = "Xval")
DF <- DF %>% pivot_longer(cols = 2:38, names_to = "Source", values_to = "Value")
DF_pvalue <- rownames_to_column(pvalue_module, var = "Xval")
DF_pvalue <- DF_pvalue %>% pivot_longer(cols = 2:38, names_to = "Source", values_to = "Value")
DF <- cbind(DF,DF_pvalue)
DF <- DF[,-c(4,5)]
names(DF) <- c("from","to","scc","pvalue")
DF$np <- ifelse(DF$scc>=0,"Positive","Negative")
DF$connect <- ifelse((abs(DF$scc)>0.2 & DF$pvalue<0.05),1,0)
links <- DF[which(DF$connect==1),]
nodes <- c(links$from,links$to)
nodes <- unique(nodes)
nodes <- data.frame(id = nodes)
nodes$degree <- 0 
test <- rbind(data.frame(table(links$from)),data.frame(table(links$to)))

matchh <-match(nodes$id,test$Var1)
nodes$degree <- test$Freq[matchh]
nodes$class <- 1
nodes$class[1:6] <-2
nodes$module_num <- 0
nodes$name <- ""
gene_m_c <- data.frame(table(moduleColors))
genus_m_c <- data.frame(table(moduleColors_g))
for (i in 1:length(nodes$module_num)) {
  if(nodes$class[i]==1){
    temp <- strsplit(nodes$id[i],"E")
    temp <- temp[[1]][2]
    nodes$module_num[i] <- gene_m_c$Freq[which(gene_m_c$moduleColors==temp)]
    nodes$name[i] <- temp
    }
  else{
    temp <- strsplit(nodes$id[i],"E")
    temp <- temp[[1]][2]
    temp <- strsplit(temp,"_")
    temp <- temp[[1]][1]
    nodes$module_num[i] <- genus_m_c$Freq[which(genus_m_c$moduleColors_g==temp)]
    nodes$name[i] <- paste(temp,"_g",sep = "")
  }
}

library(igraph)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
class(net)
net 

# We can look at the nodes, edges, and their attributes:
E(net)
V(net)

plot(net, edge.arrow.size=.4,vertex.label=paste(nodes$id,"\n(",nodes$module_num,")",sep = ""),layout=layout_components)

# Removing loops from the graph:
net <- simplify(net, remove.multiple = F, remove.loops = T) 


colrs <- c("#FF9933", "gold")
V(net)$color <- colrs[V(net)$class]
# Set node size based on audience size:
V(net)$size <- V(net)$degree*3

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label.color <- "black"

# Set edge width based on weight:
E(net)$width <- abs(E(net)$scc)*10

#change arrow size and edge color:
edge.color <- ifelse(E(net)$np=="Positive","#FF6666","#66B2FF")
pdf(file = "D:/Research/Script/WGCNA/Rplot2.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
plot(net,vertex.label.cex=.8,layout=layout_components,edge.color=edge.color,
     vertex.label=paste(nodes$name,"\n(",nodes$module_num,")",sep = "")) 
legend(x=-1.1, y=-0.8, c("Gene module","Genus module","Positive correlation","Negative correlation"), pch=c(21, 21, NA,NA),
       lty = c(NA,NA,1,1),
       col=c("#777777","#777777","#FF6666","#66B2FF"), pt.bg=colrs, pt.cex=2.5, bty="n", ncol=1)
dev.off()

plot(net,vertex.label.cex=.8,layout=layout.circle,edge.color=edge.color) 
legend(x=-1.1, y=-1.1, c("Gene module","Genus module"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2.5, bty="n", ncol=1)


ceb <- cluster_edge_betweenness(net)
dendPlot(ceb, mode="hclust")
plot(ceb, net)

####Detail module network####
# ------->> Load scc and pvalue table --------
SCC <- read.csv("D:/Research/Data/Redo/Correlation_SCC_origin.csv", header = TRUE,row.names = 1)
SCC_pvalue <- read.csv("D:/Research/Data/Redo/Correlation_pvalue_origin.csv", header = TRUE,row.names = 1)
DF <- rownames_to_column(SCC, var = "Xval")
DF <- DF %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")
DF_pvalue <- rownames_to_column(SCC_pvalue, var = "Xval")
DF_pvalue <- DF_pvalue %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")

DF <- cbind(DF,DF_pvalue)
DF <- DF[,-c(4,5)]
names(DF) <- c("from","to","scc","pvalue")
DF$np <- ifelse(DF$scc>=0,"Positive","Negative")
DF$connect <- ifelse((abs(DF$scc)>0.3 & DF$pvalue<0.05),1,0)
# ------->> genus community extraction --------
genusmodule <- cbind(moduleColors_g,moduleLabels_g)
genusmodule <- data.frame(genusmodule)
genusmodule$genus <- names(datExpr_g)

black_g <- genusmodule$genus[which(genusmodule$moduleColors_g=="black")]
blue_g <- genusmodule$genus[which(genusmodule$moduleColors_g=="blue")]
green_g <- genusmodule$genus[which(genusmodule$moduleColors_g=="green")]
red_g <- genusmodule$genus[which(genusmodule$moduleColors_g=="red")]
yellow_g <- genusmodule$genus[which(genusmodule$moduleColors_g=="yellow")]
pink_g <- genusmodule$genus[which(genusmodule$moduleColors_g=="pink")]

# ------->> gene module extraction --------
genemodule <- cbind(moduleColors,moduleLabels)
genemodule <- data.frame(genemodule)
genemodule$gene <- names(datExpr)

orange <- genemodule$gene[which(genemodule$moduleColors=="orange")]
lightyellow <- genemodule$gene[which(genemodule$moduleColors=="lightyellow")]
tan <- genemodule$gene[which(genemodule$moduleColors=="tan")]
brown4 <- genemodule$gene[which(genemodule$moduleColors=="brown4")]
saddlebrown <- genemodule$gene[which(genemodule$moduleColors=="saddlebrown")]
cyan <- genemodule$gene[which(genemodule$moduleColors=="cyan")]
lightsteelblue1 <- genemodule$gene[which(genemodule$moduleColors=="lightsteelblue1")]
violet <- genemodule$gene[which(genemodule$moduleColors=="violet")]
darkturquoise <- genemodule$gene[which(genemodule$moduleColors=="darkturquoise")]
turquoise <- genemodule$gene[which(genemodule$moduleColors=="turquoise")]
darkmagenta <- genemodule$gene[which(genemodule$moduleColors=="darkmagenta")]
brown <- genemodule$gene[which(genemodule$moduleColors=="brown")]
black <- genemodule$gene[which(genemodule$moduleColors=="black")]
floralwhite <- genemodule$gene[which(genemodule$moduleColors=="floralwhite")]
red <- genemodule$gene[which(genemodule$moduleColors=="red")]

# ------->> drawing module correlation --------
#&(abs(DF$scc)>0.35)
links2 <- DF[which((DF$connect==1) & (DF$from %in% lightyellow)&(DF$to %in%pink_g)),]
#rename genus and than can draw on the network
genusname <- links2$to
for (j in 1:length(genusname)) {
  temp <- strsplit(genusname[j],split = "\\.")
  links2$to[which(links2$to==genusname[j])] <- temp[[1]][length(temp[[1]])]
}

nodes2 <- c(links2$from,links2$to)
nodes2 <- unique(nodes2)
nodes2 <- data.frame(id = nodes2)
nodes2$degree <- 0 
test <- rbind(data.frame(table(links2$from)),data.frame(table(links2$to)))

matchh <-match(nodes2$id,test$Var1)
nodes2$degree <- test$Freq[matchh]
nodes2$class <- 1
nodes2$class[77:92] <-2

#draw network
net <- graph_from_data_frame(d=links2, vertices=nodes2, directed=F) 
class(net)
net 

# Removing loops from the graph:
net <- simplify(net, remove.multiple = F, remove.loops = T) 


colrs <- c("#FF9933", "gold")
V(net)$color <- colrs[V(net)$class]
# Set node size based on audience size:
V(net)$size <- V(net)$degree

# The labels are currently node IDs.
# Setting them to NA will render no labels:
V(net)$label.color <- "black"

# Set edge width based on weight:
E(net)$width <- abs(E(net)$scc)*5

#change arrow size and edge color:
edge.color <- ifelse(E(net)$np=="Positive","#FF6666","#66B2FF")
labell <- c()
for (i in 1:length(V(net)$degree)) {
  if(V(net)$degree[i] >=3){
    labell[i] <- nodes2$id[i]
  }else{
    labell[i] <-NA
  }
}

pdf(file = "D:/Research/Plot/1.WGCNA-BlackGenusOrangeGene1.pdf",   # The directory you want to save the file in
    width = 25, # The width of the plot in inches
    height = 25)

plot(net,vertex.label.cex=1.5,layout=layout_nicely,edge.color=edge.color,vertex.label=labell) 
dev.off()


hist(links2$scc)
net.sp <- delete_edges(net, E(net)[abs(scc)<0.31])
pdf(file = "D:/Research/Plot/1.WGCNA-BlackGenusOrangeGene1.pdf",   # The directory you want to save the file in
    width = 15, # The width of the plot in inches
    height = 15)
plot(net.sp,vertex.label.cex=.8,layout=layout_on_sphere,edge.color=edge.color,vertex.label=labell)
legend(x=-1.1, y=-1.1, c("Gene module","Genus module"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2.5, bty="n", ncol=1)
dev.off()

