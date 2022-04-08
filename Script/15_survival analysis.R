library("survival")
library("survminer")
library("readxl")

#Read data
tcga <-  read.csv("D:/Research/TCGA/TCGA.csv", header = TRUE)
tcga_cdr <-  read_excel("D:/Research/Data/Redo/TCGA-CDR-SupplementalTableS1.xlsx",sheet = "TCGA-CDR")
genus <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE)

#Grab the 138 sample from TCGA barcode and give sample name
barcode <- tcga$barcode
tcga_select <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%barcode),c("bcr_patient_barcode","OS","OS.time")]
tcga <- tcga[order(tcga$barcode),]
tcga_select<-tcga_select[order(tcga_select$bcr_patient_barcode),]
tcga_select$X <- tcga$X


#Deal with the low expression genus
genus <- genus[,-c(1286:1288)]
rownames(genus) <- genus$X
genus <- genus[,2:1285]
mean_normalized_data <- apply(genus, 2, mean)
sum(mean_normalized_data > 0)
cut_genus <- genus[, mean_normalized_data > 0]


#cbind tcga and genus data
tcga_select <- tcga_select[order(tcga_select$X),]
cut_genus <- cut_genus[order(rownames(cut_genus)),]

#rename the genus name by genus level
colname_genus <- colnames(cut_genus)
for (j in 1:length(colname_genus)) {
  if (colname_genus[j]!='X'){
    temp <- strsplit(colname_genus[j],split = "\\.")
    names(cut_genus)[names(cut_genus)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
  }
}

ana_table <- cbind(tcga_select,cut_genus)



ana_table$OS <- ifelse(ana_table$OS.time >= 2000, 0,ana_table$OS)
ana_table$OS.time <- ifelse(ana_table$OS.time >= 2000, 2000, ana_table$OS.time)
#res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
#res.cox
#summary(res.cox)

covariates <- names(cut_genus)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = ana_table)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })

res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- data.frame(res)
summ <- res[which(res$p.value<0.05),]

select_genus <- rownames(summ)
select_genus <- paste(select_genus,collapse = "+")
formular <- as.formula(paste('Surv(OS.time, OS)~', select_genus))

res.cox <- coxph(formular, data =  ana_table)
multiple_result <- summary(res.cox)

test <- survfit(res.cox)
ggsurvplot(test,data=ana_table, color = "#2E9FDF",
           ggtheme = theme_minimal())

#23 genus survival plot
res.cox <- coxph(Surv(OS.time, OS) ~ g__Ureaplasma, data = ana_table)
ggsurvplot(survfit(res.cox),data=ana_table,
           ggtheme = theme_minimal(),title = "g__Ureaplasma")


sex_df <- with(lung,
               data.frame(sex = c(1, 2), 
                          age = rep(mean(age, na.rm = TRUE), 2),
                          ph.ecog = c(1, 1)
               )
)
sex_df

fit <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit, data=lung,conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),
           ggtheme = theme_minimal())

survfit(Surv(OS.time, OS) ~ g__Ureaplasma, data = ana_table)

library(dplyr)
new <- ana_table %>% mutate(polar_exp = ifelse(g__Polaribacter >=median(ana_table$g__Polaribacter), "High", "Low"),
                            Thermo_exp = ifelse(g__Thermocrispum >=median(ana_table$g__Thermocrispum), "High", "Low"))
new$polar_exp <- factor(new$polar_exp)
new$Thermo_exp <- factor(new$Thermo_exp)

fit1 <- survfit(Surv(OS.time, OS) ~ Thermo_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE)


new1 <- new
new1$polar_exp <- "No"
new1$Thermo_exp <- "No"
for (i in 1:length(new1$g__Polaribacter)) {
  if(new1$g__Polaribacter[i]<=quantile(new1$g__Polaribacter)[2]){
    new1$polar_exp[i] <- "Low"
  }
  else if(new1$g__Polaribacter[i]>=quantile(new1$g__Polaribacter)[4]){
    new1$polar_exp[i] <-"High"
  }
}
for (i in 1:length(new1$g__Thermocrispum)) {
  if(new1$g__Thermocrispum[i]<=quantile(new1$g__Thermocrispum)[2]){
    new1$Thermo_exp[i] <- "Low"
  }
  else if(new1$g__Thermocrispum[i]>=quantile(new1$g__Thermocrispum)[4]){
    new1$Thermo_exp[i] <-"High"
  }
}
polar <- new1[-which(new1$polar_exp=="No"),]
thermo <- new1[-which(new1$Thermo_exp=="No"),]


fit1 <- survfit(Surv(OS.time, OS) ~ polar_exp, data = polar)
ggsurvplot(fit1, data = polar, pval = TRUE)


fit1 <- survfit(Surv(OS.time, OS) ~ Thermo_exp, data = thermo)
ggsurvplot(fit1, data = thermo, pval = TRUE)

#heatmap by special genus
library(pheatmap)
genus_abundance <- read.csv("D:/Research/Data/Raw/LUAD_WGS_abundance_only138.csv", header = TRUE,row.names = 1)


colname_genus <- colnames(genus_abundance)
for (j in 1:length(colname_genus)) {
  temp <- strsplit(colname_genus[j],split = "\\.")
  names(genus_abundance)[names(genus_abundance)==colname_genus[j]] <- temp[[1]][length(temp[[1]])]
}
genus <- rownames(summ)

filter <- multiple_result$coefficients[,5] <= 0.05
genus_abundance <- genus_abundance[,filter]
#names(genus_abundance)%in%genus
breaks = c(2,1.5,1,0.5,0,-0.5,-1,-1.5,-2)
test <- 200:-200
test <- test*0.01
draw <- data.frame(t(genus_abundance))
pheatmap::pheatmap(draw,scale = "row",
         color = hcl.colors(400, "RdBu"),
         breaks = test,
         show_colnames = F,
         annotation_col = subset(ana_table, select = c("OS")))




##################################################

splots <- list()
for (i in 1:length(rownames(summ))) {
  par(mar=c(10,10,3,3))
  new <- ana_table %>% mutate(exp = ifelse(ana_table[[rownames(summ)[i]]] >=median(ana_table[[rownames(summ)[i]]]), "High", "Low"))
  new$exp <- factor(new$exp)
  fit_t <- survfit(Surv(OS.time, OS) ~ exp, data = new)
  splots[[i]] <- ggsurvplot(fit_t, data = new, pval = TRUE,title = rownames(summ)[i],
             conf.int = TRUE,
             risk.table = TRUE, # Add risk table
             surv.median.line = "hv",
             risk.table.col = "exp",
             risk.table.y.text.col = T,
             palette = c("#E7B800", "#2E9FDF"),
             xlab = "Time in days",
             ggtheme = theme_light())
}
# Arrange multiple ggsurvplots and print the output
if (TRUE) {
  # Arrange and save into pdf file
  res <- arrange_ggsurvplots(splots, print = FALSE,ncol = 1, nrow = 1, risk.table.height = 0.2)
  ggsave("D:/Research/Plot/15.survivalKMplot.pdf", res, width=6, height=6)
}

summary(fit_t)$table
####store the results of multiple analysis
results <- multiple_result[["coefficients"]]
write.csv(results,file = "D:/Research/Data/Redo/Cox_surv_multipleResult.csv")

####compare with original analysis
origin_select_genus <- read.csv("D:/Research/Data/Redo/Selected_genus_in_category/prognosis_selected_genus_cutoff_2.csv",header = T)
#rename the genus name by genus level
colname_genus <- origin_select_genus$Genus
for (j in 1:length(colname_genus)) {
  if (colname_genus[j]!='X'){
    temp <- strsplit(colname_genus[j],split = "\\.")
    origin_select_genus$Genus[which(origin_select_genus$Genus==colname_genus[j])] <- temp[[1]][length(temp[[1]])]
  }
}

library(ggvenn)
x<- list(original_list = origin_select_genus$Genus,
         reanalysed_list = genus)

ggvenn(
  x, 
  fill_color = c("#FFCC99", "#99CCFF"),
  stroke_color = "#C0C0C0",
  stroke_size = 1, set_name_size = 8,
  label_sep = ",",
  text_size = 5
)
