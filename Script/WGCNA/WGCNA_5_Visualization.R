# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "D:/Research/Script/WGCNA";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Genus-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Genus-02-networkConstruction-stepByStep.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr_g, power = 3);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree_g, moduleColors_g, main = "Network heatmap plot, all genus")
plotTOM2 = 1-plotTOM
diag(plotTOM2) = 1;
TOMplot(plotTOM2, geneTree_g, moduleColors_g, main = "Network heatmap plot, all genus")

nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There¡¦s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = 0;
plotDiss = 1-plotDiss
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(MEs_g)
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengenus adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
names(MEs_g) <- paste(names(MEs_g),"_g",sep = "")

#Calculat the genus module and gene module correlation
module <- matrix(0, ncol = length(names(MEs)), nrow = length(names(MEs_g)))
module <- data.frame(module)
names(module) <- names(MEs)
rownames(module) <- names(MEs_g)
module_pvalue <- module

MEs_g$X <- rownames(MEs_g)
traitRows = match(rownames(MEs), MEs_g$X);
MEs_g = MEs_g[traitRows, -11];

MEs_g = MEs_g[order(rownames(MEs_g)),]
MEs = MEs[order(rownames(MEs)),]
for (i in 1:length(names(MEs))) {
  for (j in 1:length(names(MEs_g))) {
    m<-names(MEs)[i]
    g<-names(MEs_g)[j]
    test<-cor.test(MEs[[m]], MEs_g[[g]], method=c("spearman"))
    rho <- test$estimate[["rho"]]
    pvalue <- test$p.value
    module[which(rownames(module)==g),which(names(module)==m)] <- rho
    module_pvalue[which(rownames(module_pvalue)==g),which(names(module_pvalue)==m)] <- pvalue
  }
}
sizeGrWindow(6,10);
library(pheatmap)
pheatmap(module,main = "The correlation of Genus modules and Gene modules")

