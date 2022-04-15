setwd("D:/Research/Data/Raw")
meta <- read.csv("D:/Research/Data/Raw/Metadata_LUAD_WGS_only138.csv", header = TRUE,row.names = 1)
genus_abundance <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)

colname_genus <- colnames(genus_abundance)
for (j in 1:length(colname_genus)) {
  temp <- strsplit(colname_genus[j],split = "\\.")
  names(genus_abundance)[names(genus_abundance)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
}

target_genus <- c("g__Taupapillomavirus","g__Hypovirus","g__Betapartitivirus","g__Mamastrovirus",
                  "g__Haloplasma","g__Agrobacterium","g__Candidatus_Nitrosopelagicus","g__Polaribacter",
                  "g__Olsenella","g__Thermocrispum","g__Gemmatirosa","g__Psychrilyobacter","g__Ornithinibacillus",
                  "g__Algibacter","g__Pyramidobacter","g__Ureaplasma")

genus_matrix <- genus_abundance[,which(colnames(genus_abundance)%in%target_genus)]

#Deal with the OS time data
library("readxl")
tcga <-  read.csv("D:/Research/TCGA/TCGA.csv", header = TRUE)
tcga_cdr <-  read_excel("D:/Research/Data/Redo/TCGA-CDR-SupplementalTableS1.xlsx",sheet = "TCGA-CDR")
barcode <- tcga$barcode
tcga_select <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%barcode),c("bcr_patient_barcode","OS","OS.time")]
tcga <- tcga[order(tcga$barcode),]
tcga_select<-tcga_select[order(tcga_select$bcr_patient_barcode),]
tcga_select$X <- tcga$X

#combine data to analysis
#order the data
genus_matrix <- genus_matrix[ order(row.names(genus_matrix)), ]
meta <- meta[ order(row.names(meta)), ]
genus_matrix$stage_s <- meta$stage_severe
tcga_select <- tcga_select[order(tcga_select$X),]
genus_matrix$os.time <- tcga_select$OS.time

for (i in 1:length(genus_matrix$prognosis)) {
  if(genus_matrix$os.time[i] >= median(genus_matrix$os.time)){
    genus_matrix$prognosis[i] <- 1
  }else{
    genus_matrix$prognosis[i] <- 0
  }
}

####Logistic regression model
library(Amelia)
library(mlbench)

#separate training data and testing data
training <- genus_matrix[1:100,]
testing <- genus_matrix[101:138,]

missmap(genus_matrix, col=c("blue", "red"), legend=FALSE)

library(corrplot)
correlations <- cor(genus_matrix[,1:16])
corrplot(correlations, method="circle")

# Logistics Regression (Predict stage sever)
glm.fit <- glm(stage_s ~ g__Taupapillomavirus+g__Hypovirus+g__Betapartitivirus+g__Mamastrovirus+g__Haloplasma+
                 g__Agrobacterium+g__Candidatus_Nitrosopelagicus+g__Polaribacter+g__Olsenella+
                 g__Thermocrispum+g__Gemmatirosa+g__Psychrilyobacter+g__Ornithinibacillus+g__Algibacter+
                 g__Pyramidobacter+g__Ureaplasma, data = training, family = binomial)

summary(glm.fit)

glm.probs <- predict(glm.fit,newdata = testing,type = "response")
glm.probs[1:5]
glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
stage.test <- testing$stage_s
table(glm.pred,stage.test)
mean(glm.pred == stage.test)

library(pROC)
test_roc = roc(testing$stage_s ~ glm.probs, plot = TRUE, print.auc = TRUE)
# Logistics Regression (Predict prognosis)
glm.fit1 <- glm(prognosis ~ g__Taupapillomavirus+g__Hypovirus+g__Betapartitivirus+g__Mamastrovirus+g__Haloplasma+
                 g__Agrobacterium+g__Candidatus_Nitrosopelagicus+g__Polaribacter+g__Olsenella+
                 g__Thermocrispum+g__Gemmatirosa+g__Psychrilyobacter+g__Ornithinibacillus+g__Algibacter+
                 g__Pyramidobacter+g__Ureaplasma, data = training, family = binomial)

summary(glm.fit1)

glm.probs1 <- predict(glm.fit1,newdata = testing,type = "response")
glm.probs1[1:5]
glm.pred1 <- ifelse(glm.probs1 > 0.5, 1, 0)
prognosis.test <- testing$prognosis
table(glm.pred1,prognosis.test)
mean(glm.pred1 == prognosis.test)
test_roc = roc(testing$prognosis ~ glm.probs1, plot = TRUE, print.auc = TRUE)
roc1 <- plot.roc(testing$prognosis, glm.probs1, main="ROC comparison", percent=TRUE, col= "red", print.auc = TRUE)
roc2 <- lines.roc(testing$prognosis, glm.probs1, percent=TRUE, col="blue")
legend("bottomright", 
       legend = c("Predict1", "Predict2"), 
       col = c("red", "blue"),
       lwd = 2)
##2 genus
glm.fit <- glm(stage_s ~ g__Polaribacter+g__Thermocrispum, data = training, family = binomial)

summary(glm.fit)

glm.probs <- predict(glm.fit,newdata = training,type = "response")
glm.probs[1:5]
glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
stage.test <- training$prognosis
table(glm.pred,stage.test)
mean(glm.pred == stage.test)

#Forward selection methods(2021/12/28)
feature_list <- names(genus_matrix)
feature_list <- feature_list[1:16]
aic_list <- c()
auC_list <- c()
color_set <- c("#000000","#FF0000","#00FF00","#0000FF","#00FFFF","#FF00FF","#1E90FF","#8A2BE2",
               "#4B0082","#FFB6C1","#A0522D","#D2691E","#778899","#7FFFD4","#228B22","#FFA500")
#single feature
feature_list <- feature_list[-4]
glm.fit <- glm(prognosis ~ g__Haloplasma+g__Polaribacter+g__Algibacter, data = training, family = binomial)

glm.probs <- predict(glm.fit,newdata = testing,type = "response")
glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
prognosis.test <- testing$prognosis
test <- table(glm.pred,prognosis.test)
mean(glm.pred == prognosis.test)
roc <- plot.roc(testing$prognosis, glm.probs, main="ROC comparison", percent=TRUE, col= "#000000")#,print.auc = TRUE
auC_list[1] <- roc$auc
aic_list[1] <- glm.fit$aic
for (i in 2:length(feature_list)) {
  factor_genus <- paste(feature_list[-i], collapse="+")
  fomula <- paste("prognosis","~",factor_genus)
  glm.fit <- glm(fomula, data = training, family = binomial)
  glm.probs <- predict(glm.fit,newdata = testing,type = "response")
  glm.pred <- ifelse(glm.probs > 0.5, 1, 0)
  prognosis.test <- testing$prognosis
  roc <- lines.roc(testing$prognosis, glm.probs, percent=TRUE, col=color_set[i])
  auC_list[i] <- roc$auc
  aic_list[i] <- glm.fit$aic
}
legend("bottomright", 
       legend = feature_list, 
       col = color_set[1:13],
       lwd = 2)
auC_list <- data.frame(auC_list)
rownames(auC_list) <- feature_list
aic_list <- data.frame(aic_list)
rownames(aic_list) <- feature_list

#umap to check the clustering
library(umap)
genus_abundance <- genus_abundance[,1:1284]
umap_test <- genus_matrix[,1:16]
umap_test <- umap_test[,c(3,8,14)]
genus.umap = umap(umap_test, )

library(ggplot2)
umap_plot_df <- data.frame(genus.umap$layout)
umap_plot_df<-umap_plot_df[order(rownames(umap_plot_df)),]
umap_plot_df$stage <- meta$pathologic_stage_label
umap_plot_df$stage_severe <- meta$stage_severe
umap_plot_df$prognosis <- genus_matrix$prognosis
umap_plot_df$vital <- meta$vital_status
ggplot(
  umap_plot_df,
  aes(
    x = X1,
    y = X2,
    color = prognosis # label points with different colors for each `subgroup`
  )
) +
  geom_point()+labs(title = "Umap (Prognosis)")

genus_out <- genus_abundance
genus_out<-genus_out[order(rownames(genus_out)),]
genus_out$vital <- meta$vital_status
write.csv(genus_out,file="D:/Research/Data/Redo/genus_abundance.csv",row.names = T)

heatmap(test)
library(pheatmap)
library(grid)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.9, height=0.9, name="vp", just=c("right","top"))), action="prepend")
pheatmap(test, display_numbers = T,main = "Confusion matrix",cluster_cols = F,cluster_rows = F)
setHook("grid.newpage", NULL, "replace")
grid.text("Actual label", y=-0.07, gp=gpar(fontsize=16))
grid.text("Predict label", x=-0.07, rot=90, gp=gpar(fontsize=16))
