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

# Define numbers of genes and samples
nGenes = ncol(datExpr_g);
nSamples = nrow(datExpr_g);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_g, moduleColors_g)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits_g, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(15,15)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits_g),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               xLabelsAdj = 0.9,
               xLabelsAngle = 38,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))



#Rename genus name
genus_name <- names(datTraits)
for (j in 1:length(genus_name)) {
  temp <- strsplit(genus_name[j],split = "\\.")
  genus_name[j] <- temp[[1]][length(temp[[1]])]
}
#subgroup visualization
test_c <- moduleTraitCor[,44:63]
test_moduleTraitCor <-  moduleTraitCor[,44:63]
test_moduleTraitPvalue <-  moduleTraitPvalue[,44:63]
test_textMatrix = paste(signif(test_moduleTraitCor, 2), "\n(",
                   signif(test_moduleTraitPvalue, 1), ")", sep = "");
dim(test_textMatrix) = dim(test_moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = test_c,
               xLabels = genus_name[44:63],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = test_textMatrix,
               setStdMargins = FALSE,
               xLabelsAngle = 38,
               xLabelsAdj = 0.8,
               cex.text = 0.5,
               zlim = c(-1,1),
               font.lab.x = 0.2,
               font.lab.y = 0.2,
               main = paste("Module-trait relationships"))


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

names(datExpr)
names(datExpr)[moduleColors=="brown"]
annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")
#--------------------------------------------------------------------
##find module content
genusmodule <- cbind(moduleColors_g,moduleLabels_g)
genusmodule <- data.frame(genusmodule)
genusmodule$genus <- names(datExpr_g)
magenta <- genusmodule$genus[which(genusmodule$moduleColors_g=="magenta")]
black <- genusmodule$genus[which(genusmodule$moduleColors_g=="black")]
blue <- genusmodule$genus[which(genusmodule$moduleColors_g=="blue")]
green <- genusmodule$genus[which(genusmodule$moduleColors_g=="green")]
red <- genusmodule$genus[which(genusmodule$moduleColors_g=="red")]
yellow <- genusmodule$genus[which(genusmodule$moduleColors_g=="yellow")]
pink <- genusmodule$genus[which(genusmodule$moduleColors_g=="pink")]
brown <- genusmodule$genus[which(genusmodule$moduleColors_g=="brown")]
turquoise <- genusmodule$genus[which(genusmodule$moduleColors_g=="turquoise")]
grey <- genusmodule$genus[which(genusmodule$moduleColors_g=="grey")]

genus_G <-read.csv("D:/Research/Data/Redo/Selected_genus_in_category/gene_selected_genus_cutoff_2.csv")
genus_I <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/immune_selected_genus_cutoff_2.csv")
genus_S <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv")

tract <- match(genus_G$Genus, genusmodule$genus)
gene_module_overlap <- genusmodule[tract,]
tract <- match(genus_I$Genus, genusmodule$genus)
immune_module_overlap <- genusmodule[tract,]
tract <- match(genus_S$Genus, genusmodule$genus)
Survival_module_overlap <- genusmodule[tract,]