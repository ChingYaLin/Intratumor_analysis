# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "D:/Research/Script_Prognosis_workflow/WGCNA";
setwd(workingDir);
# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("LUAD_WGS_abundance_only138.csv",row.names = 1);
# Take a quick look at what is in the data set:
dim(femData);
names(femData);
femData <- femData[,1:1284]
datExpr0 = as.data.frame(femData);
#-------------------------------------------------------------
###Check the gene and sample witch have too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#-------------------------------------------------------------
###If gsg$allok != True, and then need to remove the sample or gene
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
#-------------------------------------------------------------
###Check outliers
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 250, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#-------------------------------------------------------------
###Loading the clinical triat data
traitData = read.csv("D:/Research/Data/Redo/Clinical Traits_AfterClean.csv",row.names = 1);
####The trait of gene####
dim(traitData)
names(traitData)
genus_G <-read.csv("D:/Research/Data/Redo/Selected_genus_in_category/gene_selected_genus_cutoff_2.csv");
genus_I <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/immune_selected_genus_cutoff_2.csv");
genus_S <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv");
genus_G <- genus_G$Genus
genus_I <- genus_I$Genus
genus_S <- genus_S$Genus
genusset <- c(genus_G,genus_I,genus_S,"X")
# remove columns that hold information we do not need.
allTraits = traitData[, which(names(traitData)%in%genusset)];
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$X);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "Gene-01-dataInput.RData")
####The trait of genus####
dim(traitData)
names(traitData)
traitData<-traitData[,-c(12)]
traitData$gender <- ifelse(traitData$gender == "FEMALE",1,0)
traitData$vital_status <- ifelse(traitData$vital_status == "Alive",1,0)
for (i in 1:length(traitData$race)) {
  if(traitData$race[i]=="WHITE"){
    traitData$race[i]<-1
  }else if(traitData$race[i]=="BLACK OR AFRICAN AMERICAN"){
    traitData$race[i]<-2
  }else{
    traitData$race[i]<-0
  }
}
traitData$pathologic_t_label <- substr(traitData$pathologic_t_label,1,1)
traitData$pathologic_n_label <- as.numeric(substr(traitData$pathologic_n_label,2,2))
for (i in 1:length(traitData$pathologic_stage_label)) {
  if(traitData$pathologic_stage_label[i]=="Stage IB" |traitData$pathologic_stage_label[i]=="Stage IA"){
    traitData$pathologic_stage_label[i]<-1
  }else if(traitData$pathologic_stage_label[i]=="Stage IIA"|traitData$pathologic_stage_label[i]=="Stage IIB"){
    traitData$pathologic_stage_label[i]<-2
  }else if(traitData$pathologic_stage_label[i]=="Stage IIIA"|traitData$pathologic_stage_label[i]=="Stage IIIB"){
    traitData$pathologic_stage_label[i]<-3}else{
    traitData$pathologic_stage_label[i]<-4
  }
}
traitData$race <- as.numeric(traitData$race)
traitData$pathologic_t_label <-as.numeric(traitData$pathologic_t_label)
traitData$pathologic_stage_label <-as.numeric(traitData$pathologic_stage_label)
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(traitData, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    main = "Sample dendrogram and trait heatmap")
datExpr_g<-datExpr
datTraits_g <- traitData
save(datExpr_g, datTraits_g, file = "Genus-01-dataInput.RData")