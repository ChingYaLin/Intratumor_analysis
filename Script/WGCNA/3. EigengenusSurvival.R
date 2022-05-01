workingDir = "D:/Research/Script/WGCNA"
setwd(workingDir)
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "Genus-01-dataInput.RData")
# Load network data saved in the second part.
lnames = load(file = "Genus-02-networkConstruction-stepByStep.RData")

nGenes = ncol(datExpr_g);
nSamples = nrow(datExpr_g);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr_g, moduleColors_g)$eigengenes
MEs = orderMEs(MEs0)
MEs = MEs[,-10]

library("survival")
library("survminer")
library("readxl")

#Read data
tcga <-  read.csv("D:/Research/TCGA/TCGA.csv", header = TRUE)
tcga_cdr <-  read_excel("D:/Research/Data/Redo/TCGA-CDR-SupplementalTableS1.xlsx",sheet = "TCGA-CDR")

#Grab the 138 sample from TCGA barcode and give sample name
barcode <- tcga$barcode
tcga_select <- tcga_cdr[which(tcga_cdr$bcr_patient_barcode%in%barcode),c("bcr_patient_barcode","OS","OS.time")]
tcga <- tcga[order(tcga$barcode),]
tcga_select<-tcga_select[order(tcga_select$bcr_patient_barcode),]
tcga_select$X <- tcga$X

#cbind tcga and eigen genus data
tcga_select <- tcga_select[order(tcga_select$X),]
MEs <- MEs[order(rownames(MEs)),]

ana_table <- cbind(tcga_select,MEs)
ana_table$OS <- ifelse(ana_table$OS.time >= 2000, 0,ana_table$OS)
ana_table$OS.time <- ifelse(ana_table$OS.time >= 2000, 2000, ana_table$OS.time)

covariates <- names(MEs)
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

select_community <- rownames(res)
select_community <- paste(select_community,collapse = "+")
formular <- as.formula(paste('Surv(OS.time, OS)~', select_community))

res.cox <- coxph(formular, data =  ana_table)
multiple_result <- summary(res.cox)

test <- survfit(res.cox)
ggsurvplot(test,data=ana_table, color = "#2E9FDF",
           ggtheme = theme_minimal())

#9 eigengenus survival plot
res.cox <- coxph(Surv(OS.time, OS) ~ MEmagenta, data = ana_table)
ggsurvplot(survfit(res.cox),data=ana_table,
           ggtheme = theme_minimal(),title = "MEmagenta")

library(dplyr)
new <- ana_table %>% mutate(magenta_exp = ifelse(MEmagenta >=median(ana_table$MEmagenta), "High", "Low"),
                            black_exp = ifelse(MEblack >=median(ana_table$MEblack), "High", "Low"),
                            brown_exp = ifelse(MEbrown >= median(ana_table$MEbrown),"High","Low"),
                            green_exp = ifelse(MEgreen >= median(ana_table$MEgreen),"High","Low"),
                            red_exp = ifelse(MEred >= median(ana_table$MEred),"High","Low"),
                            yellow_exp = ifelse(MEyellow >= median(ana_table$MEyellow),"High","Low"),
                            pink_exp = ifelse(MEpink >= median(ana_table$MEpink),"High","Low"),
                            blue_exp = ifelse(MEblue >= median(ana_table$MEblue),"High","Low"),
                            turquoise_exp = ifelse(MEturquoise >= median(ana_table$MEturquoise),"High","Low"))
new$magenta_exp <- factor(new$magenta_exp)
new$black_exp <- factor(new$black_exp)
new$brown_exp <- factor(new$brown_exp)
new$green_exp <- factor(new$green_exp)
new$red_exp <- factor(new$red_exp)
new$yellow_exp <- factor(new$yellow_exp)
new$pink_exp <- factor(new$pink_exp)
new$blue_exp <- factor(new$blue_exp)
new$turquoise_exp <- factor(new$turquoise_exp)

fit1 <- survfit(Surv(OS.time, OS) ~ magenta_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_magenta")

fit1 <- survfit(Surv(OS.time, OS) ~ black_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_black")

fit1 <- survfit(Surv(OS.time, OS) ~ brown_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_brown")

fit1 <- survfit(Surv(OS.time, OS) ~ green_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_green")

fit1 <- survfit(Surv(OS.time, OS) ~ red_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_red")

fit1 <- survfit(Surv(OS.time, OS) ~ yellow_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_yellow")

fit1 <- survfit(Surv(OS.time, OS) ~ pink_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_pink")

fit1 <- survfit(Surv(OS.time, OS) ~ blue_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_blue")

fit1 <- survfit(Surv(OS.time, OS) ~ turquoise_exp, data = new)
ggsurvplot(fit1, data = new, pval = TRUE,title = "Genus community_turquoise")
