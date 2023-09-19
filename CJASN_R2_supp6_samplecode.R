### Annotated / representative code for CJASN-2023-000878

### load and prepare CKiD data ###
### load and prepare CKiD data ###
### load and prepare CKiD data ###

library(readr)
library(dplyr)
neurocog <- read.csv("C:/Users/leeam/OneDrive - Children's Hospital of Philadelphia/Desktop/master_neurocog_long.csv")

#set covariates as factors
neurocog$sex <- as.factor(neurocog$sex)
neurocog$glom <- as.factor(neurocog$glom)
neurocog$abdelivery <- as.factor(neurocog$abdelivery)
neurocog$hypertension <- as.factor(neurocog$hypertension)
neurocog$lowses <- as.factor(neurocog$lowses)
neurocog$ancid <- as.factor(neurocog$ancid)

#log2-transform metabolite levels
neurocog[,22:643] <- log2(neurocog[,22:643])

#split based on baseline v follow-up
baseline <- subset(neurocog, visit==15) #598
followup <- subset(neurocog, visit !=15) #493

### load and prepare NiCK data ###
### load and prepare NiCK data ###
### load and prepare NiCK data ###

nick <- read_csv("nick_ckd.csv")
controls <- read_csv("nick_controls.csv")

#set covariates as factors
nick$sex <- as.factor(nick$sex)
nick$glom <- as.factor(nick$glom)
nick$abdelivery <- as.factor(nick$abdelivery)
nick$hypertension <- as.factor(nick$hypertension)


###Table -  Participant characteristics 
###Table -  Participant characteristics 
###Table -  Participant characteristics 

###baseline sample of participants with any complete neurocognitive assessments
exclude <- subset(baseline, is.na(iq))
exclude <- subset(exclude, is.na(attention))
exclude <- subset(exclude, is.na(memory))
exclude <- subset(exclude, is.na(efunction))
###filtering out participants without a complete assessment in at least 1 of 4 neurocognitive assessments

baseline <- filter(baseline, !ancid %in% exclude$ancid)

#descriptive statistics of baseline sample characteristics
summary(baseline$age)
sum(baseline$sex ==1)


###Table -  neurocognitive assessments
###Table -  neurocognitive assessments
###Table -  neurocognitive assessments

#baseline CKiD neurocognitive assessments
iq.baseline <- subset(baseline, !is.na(iq))
summary(iq.baseline$iq)
attention.baseline <- subset(baseline, !is.na(attention))
summary(attention.baseline$attention)
memory.baseline <- subset(baseline, !is.na(memory))
summary(memory.baseline$memory)
efunction.baseline <- subset(baseline, !is.na(efunction))
summary(efunction.baseline$efunction)
### repeated for CKiD follow-up, NiCK CKD, and NiCK controls


### all analyses repeated with and without adjusted for eGFR
### all analyses repeated with and without adjusted for eGFR
### all analyses repeated with and without adjusted for eGFR

### Sample linear regression analyses
### Sample linear regression analyses

#Neurocog z-score versus metabolite level
library(dplyr)
iq.baseline <- subset(baseline, !is.na(iq)) #512

estimate <- data.frame(lapply(iq.baseline[,22:643], function(x){
  summary(lm(iq~x + 
               v1b_age + sex + abdelivery +
               log2(upcr) + glom + v1b_duration + hypertension, 
             data=iq.baseline))$coef[2,1]
}))

estimate <- t(estimate)
estimate <- data.frame(estimate)

error <- data.frame(lapply(iq.baseline[,22:643], function(x){
  summary(lm(iq~x + 
               v1b_age + sex + abdelivery +
               log2(upcr) + glom + v1b_duration + hypertension, 
             data=iq.baseline))$coef[2,2]
}))

error <- t(error)
error <- data.frame(error)

p <- data.frame(lapply(iq.baseline[,22:643], function(x){
  summary(lm(iq~x + 
               v1b_age + sex + abdelivery +
               log2(upcr) + glom + v1b_duration + hypertension, 
             data=iq.baseline))$coef[2,4]
}))

p <- t(p)
p <- data.frame(p)
p$fdr <- p.adjust(p$p, method="BH")
sum(p$fdr<0.05) # FDR threshold
sum(p$p<0.05) # p<0.05

output.baseline <- data.frame(estimate, error, p)


### sample linear mixed effects analysis
### sample linear mixed effects analysis

library(lme4)
library(nlme)
library(lmerTest)

estimate <-  data.frame(lapply(iq.fu[,22:643], function(x){
  summary(lmer(iq~x + 
                 v1b_age + sex + abdelivery +
                 log2(upcr) + glom + v1b_duration + hypertension + 
                 (1|ancid), 
               data=iq.fu))$coef[2,1]
}))
estimate <- t(estimate)
estimate <- data.frame(estimate)

error <-  data.frame(lapply(iq.fu[,22:643], function(x){
  summary(lmer(iq~x + 
                 v1b_age + sex + abdelivery +
                 log2(upcr) + glom + v1b_duration + hypertension + 
                 (1|ancid), 
               data=iq.fu))$coef[2,2]
}))

error <- t(error)
error <- data.frame(error)


p <-  data.frame(lapply(iq.fu[,22:643], function(x){
  summary(lmer(iq~x + 
                 v1b_age + sex + abdelivery +
                 log2(upcr) + glom + v1b_duration + hypertension + 
                 (1|ancid), 
               data=iq.fu))$coef[2,5]
}))
p <- t(p)
p <- data.frame(p)
sum(p$p<0.05) 

output.fu <- data.frame(estimate, error, p)
output <- data.frame(output.baseline, output.fu)
write.csv(output, "output.csv")

###ratios analysis
###ratios analysis
###ratios analysis
###working memory baseline example
memory.baseline <- subset(baseline, !is.na(memory))

memory.baseline.metabolites <- select(memory.baseline, c(comp_xxxxx,
))
###select metabolites within subpathways of interest

x<-matrix(rnorm(195*47),nrow=195) #where 195=number of participant records within sample, 47=number of metabolites
x<-as.data.frame(memory.baseline.metabolites)

return<-matrix(, nrow = 195, ncol = 0)
return_name<-c()

#If we consider a/b and b/a as the same ratio, we should heve 20*19/2=190 ratios
return<-matrix(, nrow = 195, ncol = 0)
return_name<-c()
#Since we want the ratio between all paris, we should have 20*19=380 ratios
for (i in 1:dim(x)[2]){
  for (j in i:dim(x)[2]){
    if (i!=j){
      ratio=x[,i] - x[,j]
      return=cbind(return, ratio)
      return_name=c(return_name, paste(names(x)[i], "/", names(x)[j], sep=" "))
    }
  }
}
return<-as.data.frame(return)
names(return)=return_name

ratios.memory.baseline <- data.frame(memory.baseline[,1:21], return)

###calculate metabolite ratio p-value
p <-  data.frame(lapply(ratios.memory.baseline[,22:1102], function(x){
  summary(lm(memory~x + 
               v1b_age + sex + abdelivery +
               log2(upcr) + glom + v1b_duration + hypertension +log2(egfr), 
             data=ratios.memory.baseline))$coef[2,4]
}))

p <- t(p)
p <- data.frame(p)

estimate <-  data.frame(lapply(ratios.memory.baseline[,22:1102], function(x){
  summary(lm(memory~x + 
               v1b_age + sex + abdelivery +
               log2(upcr) + glom + v1b_duration + hypertension +log2(egfr), 
             data=ratios.memory.baseline))$coef[2,1]
}))

estimate <- t(estimate)
estimate <- data.frame(estimate)

output <- data.frame(estimate, p)

###calculate pgain statistic
memory.baseline <- data.frame(memory.baseline[,1:21], memory.baseline.metabolites)
#calculate individual p-values of metabolites selected for ratio testing
p <-  data.frame(lapply(memory.baseline[,22:68], function(x){
  summary(lm(memory~x + 
               age + sex + abdelivery +
               log2(upcr) + glom + duration + hypertension +log2(egfr), 
             data=memory.baseline))$coef[2,4]
}))
x<-as.data.frame(p)

return<-matrix(, nrow = 1, ncol = 0)
return_name<-c()

#calculate p-minimum for all metabolite ratio pairs
return<-matrix(, nrow = 1, ncol = 0)
return_name<-c()

for (i in 1:dim(x)[2]){
  for (j in i:dim(x)[2]){
    if (i!=j){
      minp=min(x[,i], x[,j])
      return=cbind(return, minp)
      return_name=c(return_name, paste(names(x)[i], "/", names(x)[j], sep=" "))
    }
  }
}
return<-as.data.frame(return)
names(return)=return_name
minp <- t(return)
minp <- data.frame(minp)

output <- data.frame(output, minp)
output$pgain <- output$minp / output$p
sum(output$pgain>10)
write.csv(output, "output.csv")




###Pathway enrichment analysis
phyper(#number of significant metabolites in tested subpathway,
  #number of total metabolites in tested subpathway,
  #total number of metabolites tested - number of total metabolites in tested subpathway, 
  #total number of significant metabolites, 
  lower.tail=FALSE)

