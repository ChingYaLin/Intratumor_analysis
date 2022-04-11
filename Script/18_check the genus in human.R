setwd("D:/Research/Data/Redo")
re_data <- read.csv("D:/Research/Data/Raw/GMREPO_relative_abundance_of_all_species_genus_in_all_phenotypes_summary.tsv",header = T,sep="\t")
nor_data <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)

nor_data <- nor_data[,1:1284]


min(nor_data)
mean_normalized_data <- apply(nor_data, 2, mean)
sum(mean_normalized_data > 0)
cut_genus <- nor_data[, mean_normalized_data > 0]

reference <- re_data[which(re_data$taxonomic.rank=="genus"),]
reference <-reference$scientific.name
reference <- unique(reference)



colname_genus <- colnames(cut_genus)
for (j in 1:length(colname_genus)) {
  if (colname_genus[j]!='X'){
    temp <- strsplit(colname_genus[j],split = "\\.")
    names(cut_genus)[names(cut_genus)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
  }
}
temp <- gsub("g_","",names(cut_genus))
temp <- gsub("[_]","",temp)
genus <-temp
genus <- data.frame(genus)
genus$bool <- ifelse(genus$genus%in%reference,T,F)

sum(genus$bool==T)


select_combination <- select_combination$Source
select_combination <- select_combination[!duplicated(select_combination)]

for (j in 1:length(select_combination)) {
  if (select_combination[j]!='X'){
    temp <- strsplit(select_combination[j],split = "\\.")
    select_combination[j] <- temp[[1]][length(temp[[1]])]
  }
}

select_combination <- gsub("g_","",select_combination)
select_combination <- gsub("[_]","",select_combination)

select_combination[61] <- "Ureaplasma"

select_combination <- data.frame(select_combination)
select_combination$bool <- ifelse(select_combination$select_combination%in%reference,T,F)
selectt<- select_combination[46:61,]
im_c <- c("Methylosinus","Sinorhizobium","Betabaculovirus")
im <- c(im_c , selectt$select_combination)
genus_focus <- select_combination[which(select_combination$select_combination%in%im),]

sum(select_combination$bool==T)
select_combination$select_combination[53] <- "Polaribacter"
