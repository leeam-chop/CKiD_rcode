#SAMPLE CODE

#Load Data
library(readr)
longitudinal <- read.csv("longitudinal.csv")

#variable prep
library(dplyr)
longitudinal <- data.frame(longitudinal)
longitudinal$sex <- as.factor(longitudinal$sex)
longitudinal$glom <- as.factor(longitudinal$glom)
longitudinal[,21:642] <- log(longitudinal[,21:642]) #natural log transform metabolite levels
longitudinal$timepoint <- as.factor(longitudinal$timepoint)

###Table 1 descriptive statistics
baseline <- subset(longitudinal, visit==15)
summary(baseline$v1b_age)
summary(baseline$bmiz)
summary(baseline$egfr)
summary(baseline$upcr)
summary(baseline$duration)

followup <- subset(longitudinal, visit != 15)

followed <- filter(baseline, ancid %in% followup$ancid)
summary(followed$v1b_age)
summary(followed$bmiz)
summary(followed$egfr)
summary(followed$upcr)
summary(followed$duration)

### Question 1 - What metabolites are related to eGFR in repeated cross-sectional modeling
library(lme4)
library(nlme)
library(lmerTest)

lme.egfr <-  data.frame(lapply(longitudinal[,21:642], function(x){
  summary(lmer(x ~ log(egfr) + age + sex + (1|ancid), 
               data=longitudinal))$coef[2,c(1,2,5)]
}))
lme.egfr <- t(lme.egfr)
lme.egfr <- data.frame(lme.egfr)
lme.egfr <- lme.egfr %>%
  rename("p" = "Pr...t..",
         "error" = "Std..Error")
lme.egfr$fdr <- p.adjust(lme.egfr$p, method="BH")
sum(lme.egfr$p<0.05) 
sum(lme.egfr$fdr<0.05) 

###Per 10% unit transformation
lme.egfr$E10 <- ((1.1^lme.egfr$Estimate) - 1 )*100
lme.egfr$lower95 <- lme.egfr$Estimate - 1.96*lme.egfr$error
lme.egfr$upper95 <- lme.egfr$Estimate + 1.96*lme.egfr$error
lme.egfr$E10_lower95 <- ((1.1^lme.egfr$lower95) - 1 )*100
lme.egfr$E10_upper95 <- ((1.1^lme.egfr$upper95) - 1 )*100

spearman.egfr <- apply(longitudinal[,21:642], 2, function(x) cor.test(log(longitudinal$egfr), x, method="spearman")$estimate)
spearman.egfr <- data.frame(spearman.egfr)

output <- data.frame(lme.egfr, spearman.egfr)

### Question 2a - What metabolites are related to UPCR in repeated cross-sectional modeling
library(lme4)
library(nlme)
library(lmerTest)

lme.upcr <-  data.frame(lapply(longitudinal[,21:642], function(x){
  summary(lmer(x ~ log(upcr) +age + sex + (1|ancid), 
               data=longitudinal))$coef[2,c(1,2,5)]
}))
lme.upcr <- t(lme.upcr)
lme.upcr <- data.frame(lme.upcr)
lme.upcr <- lme.upcr %>%
  rename("p" = "Pr...t..",
         "error" = "Std..Error")
lme.upcr$fdr <- p.adjust(lme.upcr$p, method="BH")
sum(lme.upcr$p<0.05) #414 metabolites
sum(lme.upcr$fdr<0.05) #405 metabolites

###Per 10% unit transformation
lme.upcr$E10 <- ((1.1^lme.upcr$Estimate) - 1 )*100
lme.upcr$lower95 <- lme.upcr$Estimate - 1.96*lme.upcr$error
lme.upcr$upper95 <- lme.upcr$Estimate + 1.96*lme.upcr$error
lme.upcr$E10_lower95 <- ((1.1^lme.upcr$lower95) - 1 )*100
lme.upcr$E10_upper95 <- ((1.1^lme.upcr$upper95) - 1 )*100

spearman.upcr <- apply(longitudinal[,21:642], 2, function(x) cor.test(log(longitudinal$upcr), x, method="spearman")$estimate)
spearman.upcr <- data.frame(spearman.upcr)

output <- data.frame(lme.upcr, spearman.upcr)

###question 3 - Do metabolites change over time between V3->V5, n=530 person visits
followup <- subset(longitudinal, visit!=15)

e1 <-  data.frame(lapply(followup[,21:642], function(x){
  summary(lmer(x ~ timepoint + v1b_age + sex + (1|ancid), 
               data=followup))$coef[2,1]
}))
e1 <- t(e1)
e1 <- data.frame(e1)

p1 <-  data.frame(lapply(followup[,21:642], function(x){
  summary(lmer(x ~ timepoint + v1b_age + sex + (1|ancid), 
               data=followup))$coef[2,5]
}))
p1 <- t(p1)
p1 <- data.frame(p1)
p1$fdr <- p.adjust(p1$p1, method="BH")


output <- data.frame(e1, p1)
sum(output$fdr<0.05) ###158/622 metabolites demonstrate significant time differences between V3->V5
sum(output$p<0.05)

####Delta-delta analysis
### How does change in metabolite level associate with change in eGFR and proteinuria in the follow-up visits
longitudinal <- subset(longitudinal, !is.na(upcr))
v3 <- subset(longitudinal, visit==30)
v5 <- subset(longitudinal, visit==50)
v3 <- filter(v3, ancid %in% v5$ancid)
v5 <- filter(v5, ancid %in% v3$ancid)

delta_metabolites <- v5[,21:642] - v3[,21:642] #% change in metabolite level
delta_egfr <- log(v5$egfr) - log(v3$egfr) #%change in eGFR level
base_egfr <- log(v3$egfr)
delta_upcr <- log(v5$upcr) - log(v3$upcr)
base_upcr <- log(v3$upcr)
delta <- data.frame(v5[,1:20], base_egfr, base_upcr, delta_egfr, delta_upcr, delta_metabolites)

summary(delta$delta_egfr)
summary(delta$delta_upcr)

output <-  data.frame(lapply(delta[,25:646], function(x){
  summary(lm(x ~ delta_egfr + delta_upcr +
               v1b_age + sex + base_egfr + base_upcr, 
             data=delta))$coef[2,c(1,2,4)]
}))
output <- t(output)
output <- data.frame(output) %>%
  rename("p" = "Pr...t..",
         "error" = "Std..Error")
output$fdr <- p.adjust(output$p, method="BH")
sum(output$p<0.05) #133 metabolites
sum(output$fdr<0.05) #52 metabolites

output$E10 <- ((1.1^output$Estimate) - 1 )*100
output$lower95 <- output$Estimate - 1.96*output$error
output$upper95 <- output$Estimate + 1.96*output$error
output$E10_lower95 <- ((1.1^output$lower95) - 1 )*100
output$E10_upper95 <- ((1.1^output$upper95) - 1 )*100

write.csv(output, "output.csv")


spearman.egfr <- apply(delta[,25:646], 2, function(x) cor.test(delta$delta_egfr, x, method="spearman")$estimate)
spearman.egfr <- data.frame(spearman.egfr)
write.csv(spearman.egfr, "output.csv")


output <-  data.frame(lapply(delta[,25:646], function(x){
  summary(lm(x ~ delta_upcr + delta_upcr +
               v1b_age + sex + base_upcr + base_upcr, 
             data=delta))$coef[2,c(1,2,4)]
}))
output <- t(output)
output <- data.frame(output) %>%
  rename("p" = "Pr...t..",
         "error" = "Std..Error")
output$fdr <- p.adjust(output$p, method="BH")
sum(output$p<0.05) #133 metabolites
sum(output$fdr<0.05) #52 metabolites

output$E10 <- ((1.1^output$Estimate) - 1 )*100
output$lower95 <- output$Estimate - 1.96*output$error
output$upper95 <- output$Estimate + 1.96*output$error
output$E10_lower95 <- ((1.1^output$lower95) - 1 )*100
output$E10_upper95 <- ((1.1^output$upper95) - 1 )*100


spearman.upcr <- apply(delta[,25:646], 2, function(x) cor.test(delta$delta_upcr, x, method="spearman")$estimate)
spearman.upcr <- data.frame(spearman.upcr)
output <- data.frame(output, spearman.upcr)
write.csv(spearman.upcr, "output.csv")

write.csv(output, "output.csv")

### Does baseline metabolite levels predict change in proteinuria over time?
v5 <- longitudinal %>%
  subset(visit ==50) 
v1b <- longitudinal %>%
  subset(visit==15)
v5 <- v5 %>% subset(ancid %in% v1b$ancid)
v1b <- v1b %>% subset(ancid %in% v5$ancid)
#260 participants with upcr and metabolomics at baseline and year 4
v1b$upcr.v5 <- v5$upcr
v1b$delta.upcr <- log(v1b$upcr.v5) - log(v1b$upcr)

output <-  data.frame(lapply(v1b[,21:642], function(x){
  summary(lm(delta.upcr ~ x + 
               v1b_age + sex + log(egfr) + log(upcr), 
             data=v1b))$coef[2,c(1,2,4)]
}))
output <- t(output)
output <- data.frame(output) %>%
  rename("p" = "Pr...t..",
         "error" = "Std..Error")
output$fdr <- p.adjust(output$p, method="BH")
sum(output$p<0.05) #78 metabolites
sum(output$fdr<0.05) #1 metabolites

output$E10 <- ((1.1^output$Estimate) - 1 )*100
output$lower95 <- output$Estimate - 1.96*output$error
output$upper95 <- output$Estimate + 1.96*output$error
output$E10_lower95 <- ((1.1^output$lower95) - 1 )*100
output$E10_upper95 <- ((1.1^output$upper95) - 1 )*100


spearman.upcr <- apply(v1b[,21:642], 2, function(x) cor.test(v1b$delta.upcr, x, method="spearman")$estimate)
spearman.upcr <- data.frame(spearman.upcr)
output <- data.frame(output, spearman.upcr)
write.csv(output, "output.csv")

