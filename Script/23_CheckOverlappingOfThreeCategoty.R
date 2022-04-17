setwd('D:/Research/Data/Redo/Selected_genus_in_category')
#Read three categories selected genus and it's prevalence information (here the prevalence is calculate by value of 2)
gene_select_2 <- read.csv("gene_selected_genus_cutoff_2.csv")
immune_select_2 <- read.csv("immune_selected_genus_cutoff_2.csv")
prognosis_select_2 <- read.csv("prognosis_selected_genus_cutoff_2.csv")
#Read three categories selected genus and it's prevalence information (here the prevalence is calculate by value of 1)
gene_select_1 <- read.csv("gene_selected_genus_cutoff_1.csv")
immune_select_1 <- read.csv("immune_selected_genus_cutoff_1.csv")
prognosis_select_1 <- read.csv("prognosis_selected_genus_cutoff_1.csv")

#Find intersect of three category's genus
shared_genus <- intersect(immune_select_2$Genus,gene_select_2$Genus)
