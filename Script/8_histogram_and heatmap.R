library(tibble)
library(ggplot2)
library(tidyr)
library(pheatmap)

####import data####
#data is all gene correlation
data <- read.csv("D:/Research/Data/Redo/Correlation_SCC_origin.csv", header = TRUE,row.names = 1)
ran_data <- read.csv("D:/Research/Data/Redo/Correlation_SCC_random.csv", header = TRUE,row.names = 1)

#construct the value to one column that can plot by histogram
DF <- rownames_to_column(data, var = "Xval")
DF <- DF %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")
DF_ran <- rownames_to_column(ran_data, var = "Xval")
DF_ran <- DF_ran %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")


#plot the histogram
ggplot(DF)+
  geom_histogram(aes(x=Value),binwidth = 0.05,color = "royalblue",fill="slategray1")+
  labs(title = "Correlation distribution",subtitle = "15323 genes with 1146 genus")+
  xlab("Spearman Correlation coefficient")
hist(DF$Value,main = "Correlation distribution",
     col = "slategray1",
     border = "royalblue",
     xlab = "Spearman Correlation coefficient")

ggplot(DF_ran)+
  geom_histogram(aes(x=Value),binwidth = 0.05,color = "seagreen",fill="lightgreen")+
  labs(title = "Correlation distribution (random permutation)",subtitle = "15323 genes with 1146 genus")+
  xlab("Spearman Correlation coefficient")


####plot the heatmap(correlation)####
#rename the genus first
colname_wgs <- colnames(data)
for (j in 1:length(colname_wgs)) {
  temp <- strsplit(colname_wgs[j],split = "\\.")
  names(data)[names(data)==colname_wgs[j]] <- temp[[1]][length(temp[[1]])]
}
heatmap(t(data),main = 'Heatmap(origin)')

colname_wgs <- colnames(ran_data)
for (j in 1:length(colname_wgs)) {
  temp <- strsplit(colname_wgs[j],split = "\\.")
  names(ran_data)[names(ran_data)==colname_wgs[j]] <- temp[[1]][length(temp[[1]])]
}
heatmap(t(ran_data),main = 'Heatmap(random)')

res <- wilcox.test(DF$Value, DF_ran$Value, alternative = "two.sided")
boxplot(DF$Value,DF_ran$Value)


#plot density
DF$class <- "Origin"
DF_ran$class <- "Random permutation"
test <- rbind(DF,DF_ran)
pdf("D:/Research/Plot/8.SccDensity_WGS.pdf",
    width = 10,
    height = 8)
par(mar=c(3,3,3,3))
ggplot(test)+
  geom_density(aes(x=Value,color = class,fill = class),alpha = 0.2)+
  labs(title = "SCC Density")
dev.off()

count_origin <- table(cut(DF$Value,breaks = c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)))
barplot(count_origin)
count_random <- table(cut(DF_ran$Value,breaks = c(-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)))
barplot(count_random)

#plot the difference count between the origin and random group
count_dif <- count_origin-count_random
barplot(count_dif,col = "tan1",border="darkorange",main="The barplot of origin minus random",xlab="SCC interval",ylab="count")

summary_count <- data.frame(cbind(count_origin,count_random,count_dif))



#select less gene and genus to plot density
fil_origin <- DF[which((DF$Xval %in% cut_gene)&(DF$Source %in% cut_genus)),]
fil_random <- DF_ran[which((DF_ran$Xval %in% cut_gene)&(DF_ran$Source %in% cut_genus)),]
test <- rbind(fil_origin,fil_random)
ggplot(test)+
  geom_density(aes(x=Value,color = class,fill = class),alpha = 0.2)+
  labs(title = "Density")

######Immune
data <- read.csv("D:/Research/Data/Redo/Immune_Correlation_SCC_origin.csv", header = TRUE,row.names = 1)
ran_data <- read.csv("D:/Research/Data/Redo/Immune_Correlation_SCC_random.csv", header = TRUE,row.names = 1)
DF <- rownames_to_column(data, var = "Xval")
DF <- DF %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")
DF_ran <- rownames_to_column(ran_data, var = "Xval")
DF_ran <- DF_ran %>% pivot_longer(cols = 2:1147, names_to = "Source", values_to = "Value")

hist(DF$Value,main = "Correlation distribution",
     col = "slategray1",
     border = "royalblue",
     xlab = "Spearman Correlation coefficient")


select_combination <- DF[which(abs(DF$Value)>=0.3),]
select_combination <- select_combination[order(select_combination$Source),]

genus <- select_combination$Source
genus_rm <- genus[!duplicated(genus)]
test <- table(genus)
test <- data.frame(test)


#### Scatter plot
filter_genus <- test[which(test$Freq>=3),]
filter <- select_combination[which(select_combination$Source%in%filter_genus$genus),]
filter <- filter[which(abs(filter$Value)>=0.3),]
immune_data <- read.csv("D:/Research/Data/Redo/estimation_matrix.csv", header = TRUE,row.names = 1)
mean_check <- apply(immune_data, 1, mean)
hist(mean_check)
sum(mean_check > 0)
immune_data <- immune_data[mean_check > 0, ]

genus_abundance <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)
genus_abundance <- genus_abundance[,which(names(genus_abundance)%in%filter_genus$genus)]
genus_abundance <- data.frame(t(genus_abundance))

rownames(immune_data) <- rownames(data)

immune_data <- immune_data[which(rownames(immune_data)%in%filter$Xval),]
immune_data <- immune_data[,order(names(immune_data))]
genus_abundance <- genus_abundance[,order(names(genus_abundance))]

plot_data <- rbind(immune_data,genus_abundance)
plot_data <- data.frame(t(plot_data))
library("ggpubr")
ggscatter(plot_data, x = "T.cell.NK_XCELL", y = "k__Viruses.f__Baculoviridae.g__Betabaculovirus", 
          add = "loess", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab ="T.cell.NK_XCELL", ylab = "g__Betabaculovirus")+
  labs(title="Spearman correlation")