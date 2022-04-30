workingDir = "D:/Research/Script/WGCNA";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Genus-01-dataInput.RData");
# Load network data saved in the second part.
lnames = load(file = "Genus-02-networkConstruction-stepByStep.RData");
lnames
immune <- read.csv("estimation_matrix.csv",header = T,row.names = 1)
immune <- data.frame(t(immune))
cibersort_abs <- immune[,29:50]
cibersort_abs <- cibersort_abs[order(rownames(cibersort_abs)),]

nGenes = ncol(datExpr_g);
nSamples = nrow(datExpr_g);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_g, moduleColors_g)$eigengenes
MEs = orderMEs(MEs0)
MEs = MEs[,-10]
cibersort_abs <- cibersort_abs[,-5]
#rename column name
immunename <- names(cibersort_abs)
for (j in 1:length(immunename)) {
  temp <- strsplit(immunename[j],split = "_")
  names(cibersort_abs)[j] <- temp[[1]][1]
}

moduleimmuneCor = cor(MEs, cibersort_abs, use = "p",method = "spearman");
moduleimmunePvalue = corPvalueStudent(moduleimmuneCor, nSamples);


sizeGrWindow(25,15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleimmuneCor, 2), "\n(",
                   signif(moduleimmunePvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleimmuneCor)

pdf(file = "D:/Research/Script/WGCNA/Rplot1.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6)
par(mar = c(10, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleimmuneCor,
               xLabels = names(cibersort_abs),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               xLabelsAngle = 40,
               zlim = c(-1,1),
               main = paste("Community-immune infiltration relationships"))
dev.off()

hist(moduleimmuneCor,xlab = "SCC",main = "Histogram of Community Correlation")


####random####
ran_cibersort <- cibersort_abs
for (i in 1:length(names(ran_cibersort))) {
  ran_cibersort[[names(ran_cibersort)[i]]] <- sample(ran_cibersort[[names(ran_cibersort)[i]]])
  
  
}


moduleimmuneCor_r = cor(MEs, ran_cibersort, use = "p",method = "spearman");
moduleimmunePvalue_r = corPvalueStudent(moduleimmuneCor_r, nSamples);


sizeGrWindow(25,15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleimmuneCor_r, 2), "\n(",
                   signif(moduleimmunePvalue_r, 1), ")", sep = "");
dim(textMatrix) = dim(moduleimmuneCor_r)

par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleimmuneCor_r,
               xLabels = names(ran_cibersort),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               xLabelsAdj = 0.75,
               xLabelsAngle = 30,
               zlim = c(-1,1),
               main = paste("Community-immune infiltration relationships"))

hist(moduleimmuneCor_r,xlab = "SCC",main = "Histogram of Community Correlation")


library(tibble)
library(ggplot2)
library(tidyr)
library(pheatmap)
DF <- rownames_to_column(data.frame(moduleimmuneCor), var = "Xval")
DF <- DF %>% pivot_longer(cols = 2:22, names_to = "Source", values_to = "Value")
DF_ran <- rownames_to_column(data.frame(moduleimmuneCor_r), var = "Xval")
DF_ran <- DF_ran %>% pivot_longer(cols = 2:22, names_to = "Source", values_to = "Value")

DF$class <- "Origin"
DF_ran$class <- "Random permutation"
test <- rbind(DF,DF_ran)
ggplot(test)+
  geom_density(aes(x=Value,color = class,fill = class),alpha = 0.2)+
  labs(title = "Density of correlation values")