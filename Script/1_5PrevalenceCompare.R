#### Functions ####
RenameGenus <- function(df) {
  #This function can rename the data frame that the column name is genus
  bacteria_name <- colnames(df)
  for (i in 1:length(bacteria_name)) {
    if(bacteria_name[i]!="X"){
      temp <- strsplit(bacteria_name[i],split = "\\.")
      names(df)[names(df)==bacteria_name[i]] <- temp[[1]][length(temp[[1]])]
    }
  }
  df <- df[,!(names(df) %in% c("X"))]
  return(df)
}

Prevalence_cal <- function(df,cutoff=0){
  # This function can calculate the genus prevalence of whole data frame
  ## Input: df--data frame that column is genus, raw is sample,
  ##        cutoff--the baseline of genus that seen presented in the sample, default is 0
  ## Output: a data frame that have whole genus prevalence
  prevalence <- colSums(df >= cutoff)
  prevalence <- prevalence/length(df[,1])
  return(data.frame(prevalence))
}

####Main####
# ----> Setting environment and load files----
setwd("D:/Research/Data/Raw")
### Whole LUAD sample compare(WGS:162 samples,RNA:519 samples)
#df_wgs <- read.csv("LUAD_WGS_abundance.csv") 
#df_rna <- read.csv("LUAD_RNA_seq_abundance.csv")
### Same LUAD sample compare(WGS:138 samples,RNA:138 samples)
df_wgs <- read.csv("LUAD_WGS_abundance_only138.csv")
df_rna <- read.csv("LUAD_RNA_seq_abundance_only138.csv")

# ----> Rename the column's name to genus level----
df_wgs <- RenameGenus(df_wgs)
df_rna <- RenameGenus(df_rna)

# ----> Remove the column of contaminant----
contaminant <- c("contaminant1Harvard","contaminant2HarvardCanadaBaylorWashU","contaminant3AllSeqCenters")
df_rna = df_rna[,!(names(df_rna) %in% contaminant)]
df_wgs = df_wgs[,!(names(df_wgs) %in% contaminant)]

# ----> Calculate the prevalence----
prevalence_rna <- Prevalence_cal(df_rna,2)
prevalence_wgs <- Prevalence_cal(df_wgs,2)

# ----> Combine result----
compare_table <- data.frame(prevalence_rna,prevalence_wgs)
names(compare_table) <- c("RNA","WGS")

# ----> Visualization----
library(ggplot2)
# Sort the value
compare_table <- compare_table[order(-compare_table$RNA,-compare_table$WGS),]
num <- c(1:length(compare_table$RNA))
compare_table$num <- num
compare_table$compare <- compare_table$RNA-compare_table$WGS
# Draw the line
ggplot(compare_table) +
  geom_histogram(aes(compare),fill = "blue",binwidth = 0.01)+
  labs(title = "Prevalence compare", subtitle = "RNA-seq prevalence - WGS prevalence",x="value", y="count")
new <- compare_table[,c("num","compare")]
new <- new[order(new$compare),]
new$num <- num
ggplot(new,aes(x=num,y=compare))+
  geom_line(color="blue")+
  geom_hline(yintercept = 0,color = "red")+
  labs(title = "Prevalence compare", subtitle = "RNA-seq prevalence - WGS prevalence",x="microbiome", y="comparison")


#count the number of > and <
rna_bigger <- sum(compare_table$compare>=0.1)
wgs_bigger <- sum(compare_table$compare<=(-0.1))
equal <- sum(compare_table$compare>(-0.1) & compare_table$compare <0.1)
cat("RNA>0.1:",rna_bigger,"WGS>0.1:",wgs_bigger,"Equal:",equal)
