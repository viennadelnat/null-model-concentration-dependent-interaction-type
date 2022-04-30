#Whether warming magnifies the toxicity of a pesticide is strongly dependent on the concentration and the null model
#Vienna Delnat, Lizanne Janssens and Robby Stoks
#Aquatic Toxicology (2019)
#R code tested on 23/03/2022


#####Packages#####

install.packages("car")     
install.packages("lme4")     
install.packages("lsmeans")     
install.packages("effects")     
install.packages("afex")

library(car)     
library(lme4)
library(lsmeans)
library(effects)     
library(afex)

##Check versions of packages and R
sessionInfo()

##Set working directory to source file
#RStudio -> Session -> Set Working Directory...-> To Source File Location


#####Datasets#####

#Dataset of survival
dataBIN=read.csv("Delnat-et-al_SurvivalBinomial.csv", sep=",", na.strings=c(""))
#Set correct data types 
Factors <- c("Temperature", "CPF", "Replicate", "EggDate", "StartDate")
dataBIN[Factors] <- do.call(cbind.data.frame, lapply(dataBIN[Factors], as.factor))
str(dataBIN)

#Dataset of growth rate
data=read.csv("Delnat-et-al_GrowthRate.csv", sep=",", na.strings=c(""))
#Set correct data types 
Factors <- c("Temperature", "CPF", "Replicate", "EggDate", "StartDate")
data[Factors] <- do.call(cbind.data.frame, lapply(data[Factors], as.factor))
str(data)

#Dataset of physiological endpoints
dataPhys=read.csv("Delnat-et-al_Physiology.csv", sep=",", na.strings=c(""))
#Set correct data types 
Factors <- c("Temperature", "CPF", "Jar")
dataPhys[Factors] <- do.call(cbind.data.frame, lapply(dataPhys[Factors], as.factor))
str(dataPhys)


#####Survival#####

#Effect coding due to interaction effects (also use Anova type III instead of type II)
set_sum_contrasts()

#Generalized linear mixed model with a binomial error structure and the logit link
#To take into account that animals from the same vial were not independent we added vial (=Replicate) to the model as a random factor.
glmerSurvivalEC=glmer(SurvivalEC ~ Temperature*CPF + (1|Replicate), data=dataBIN, na.action=na.omit, family=binomial(link=logit),
                 control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerSurvivalEC, type="III") 

#Posthoc - Contrast analysis
interact1<-pairs(lsmeans(glmerSurvivalEC, ~Temperature|CPF, adjust="none"))
interact2<-pairs(lsmeans(glmerSurvivalEC, ~CPF|Temperature, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#lsmeans and standard errors used to make the figure in Excel
lsmeans(glmerSurvivalEC, ~Temperature*CPF, adjust="none", type="response")

#Assumption - Dispersion parameter
glmSurvivalEC=glm(SurvivalEC ~ Temperature*CPF, data=dataBIN, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmSurvivalEC) 


#####Growth rate#####

#general linear model with a normal error structure and the identity link
lmGrowth=lm(GrowthDev ~ Temperature*CPF, data=data, na.action=na.omit) 
Anova(lmGrowth, type="III") 

#Posthoc - Contrast analysis
pairs(lsmeans(lmGrowth, ~CPF, adjust="fdr"))

#lsmeans and standard errors used to make the figure in Excel
lsmeans(lmGrowth, ~Temperature*CPF, adjust="none")

#Outliers and influential observations
outlierTest(lmGrowth) 
cd=cooks.distance(lmGrowth)
which(cd>1) 

#Assumptions
#Normality of residuals
shapiro.test(resid(lmGrowth))                  
hist(resid(lmGrowth))   
qqnorm(resid(lmGrowth))    
qqline(resid(lmGrowth))     
#Homogeneity of variance
leveneTest(GrowthDev ~ Temperature*CPF, data = data)
aggregate(GrowthDev ~ Temperature*CPF, data = data, var)


#####Acetylcholinesterase activity#####

#general linear model with a normal error structure and the identity link
lmAChE=lm(AChE_Protein ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
dataPhys$AChE_ProteinTR=log(dataPhys$AChE_Protein)
lmAChE=lm(AChE_ProteinTR ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
Anova(lmAChE, type="III") 

#Posthoc - Contrast analysis
pairs(lsmeans(lmAChE, ~CPF, adjust="fdr"))

#lsmeans and standard errors used to make the figure in Excel
lsmeans(lmAChE, ~Temperature*CPF, adjust="none")

#Outliers and influential observations
outlierTest(lmAChE) 
cd=cooks.distance(lmAChE)
which(cd>1) 

#Assumptions
#Normality of residuals
shapiro.test(resid(lmAChE))                  
hist(resid(lmAChE))    
qqnorm(resid(lmAChE))    
qqline(resid(lmAChE))     
#Homogeneity of variance
leveneTest(AChE_ProteinTR ~ Temperature*CPF, data = dataPhys)
aggregate(AChE_ProteinTR ~ Temperature*CPF, data = dataPhys, var)
leveneTest(AChE_ProteinTR ~ Temperature*CPF, data = dataPhys)
aggregate(AChE_ProteinTR ~ Temperature*CPF, data = dataPhys, var)


#####Total fat content#####

#general linear model with a normal error structure and the identity link
lmTotalFat=lm(microgTotalFat_mgL4 ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
dataPhys$microgTotalFat_L4TR=log(dataPhys$microgTotalFat_mgL4)
lmTotalFat=lm(microgTotalFat_L4TR ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
Anova(lmTotalFat, type="III") 

#Posthoc - Contrast analysis
pairs(lsmeans(lmTotalFat, ~CPF, adjust="none"))

#lsmeans and standard errors used to make the figure in Excel
lsmeans(lmTotalFat, ~Temperature*CPF, adjust="none")

#Outliers and influential observations
outlierTest(lmTotalFat) 
cd=cooks.distance(lmTotalFat)
which(cd>1) 

#Assumptions
#Normality of residuals
shapiro.test(resid(lmTotalFat))                  
hist(resid(lmTotalFat))    
qqnorm(resid(lmTotalFat))    
qqline(resid(lmTotalFat))     
#Homogeneity of variance
leveneTest(microgTotalFat_mgL4 ~ Temperature*CPF, data = dataPhys)
aggregate(microgTotalFat_mgL4 ~ Temperature*CPF, data = dataPhys, var)
leveneTest(microgTotalFat_L4TR ~ Temperature*CPF, data = dataPhys)
aggregate(microgTotalFat_L4TR ~ Temperature*CPF, data = dataPhys, var)


#####Electron transport system activity#####

#general linear model with a normal error structure and the identity link
lmOxygenConsumption=lm(O2_protein ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
Anova(lmOxygenConsumption, type="III") 

#Posthoc - Contrast analysis
interact1<-pairs(lsmeans(lmOxygenConsumption, ~Temperature|CPF, adjust="none"))
interact2<-pairs(lsmeans(lmOxygenConsumption, ~CPF|Temperature, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#lsmeans and standard errors used to make the figure in Excel
lsmeans(lmOxygenConsumption, ~Temperature*CPF, adjust="none")

#Outliers and influential observations
outlierTest(lmOxygenConsumption) 
cd=cooks.distance(lmOxygenConsumption)
which(cd>1) 

#Assumptions
#Normality of residuals
shapiro.test(resid(lmOxygenConsumption))                  
hist(resid(lmOxygenConsumption))    
qqnorm(resid(lmOxygenConsumption))    
qqline(resid(lmOxygenConsumption))     
#Homogeneity of variance
leveneTest(O2_protein ~ Temperature*CPF, data = dataPhys)
aggregate(O2_protein ~ Temperature*CPF, data = dataPhys, var)


#####Superoxide anion concentration#####

#general linear model with a normal error structure and the identity link
lmSuperoxideAnion=lm(SuperoxideAnion_mgL4 ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
dataPhys$SuperoxideAnion_mgL4TR=log(dataPhys$SuperoxideAnion_mgL4)
lmSuperoxideAnion=lm(SuperoxideAnion_mgL4TR ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
Anova(lmSuperoxideAnion, type="III") 

#lsmeans and standard errors used to make the figure in Excel
lsmeans(lmSuperoxideAnion, ~Temperature*CPF, adjust="none")

#Outliers and influential observations
outlierTest(lmSuperoxideAnion) 
cd=cooks.distance(lmSuperoxideAnion)
which(cd>1) 

#Assumptions
#Normality of residuals
shapiro.test(resid(lmSuperoxideAnion))                  
hist(resid(lmSuperoxideAnion))    
qqnorm(resid(lmSuperoxideAnion))    
qqline(resid(lmSuperoxideAnion))     
#Homogeneity of variance
leveneTest(SuperoxideAnion_mgL4 ~ Temperature*CPF, data = dataPhys)
aggregate(SuperoxideAnion_mgL4 ~ Temperature*CPF, data = dataPhys, var)
leveneTest(SuperoxideAnion_mgL4TR ~ Temperature*CPF, data = dataPhys)
aggregate(SuperoxideAnion_mgL4TR ~ Temperature*CPF, data = dataPhys, var)


#####Malondialdehyde level#####

#general linear model with a normal error structure and the identity link
lmMDA=lm(MDA_Fat ~ Temperature*CPF, data=dataPhys, na.action=na.omit) 
Anova(lmMDA, type="III") 

#Posthoc - Contrast analysis
pairs(lsmeans(lmMDA, ~CPF, adjust="fdr"))

#lsmeans and standard errors used to make the figure in Excel
lsmeans(lmMDA, ~Temperature*CPF, adjust="none")

#Outliers and influential observations
outlierTest(lmMDA) 
cd=cooks.distance(lmMDA)
which(cd>1) 

#Assumptions
#Normality of residuals
shapiro.test(resid(lmMDA))                  
hist(resid(lmMDA))    
qqnorm(resid(lmMDA))    
qqline(resid(lmMDA))     
#Homogeneity of variance
leveneTest(MDA_Fat ~ Temperature*CPF, data = dataPhys)
aggregate(MDA_Fat ~ Temperature*CPF, data = dataPhys, var)


######Save Rdata######
save.image(file="Rdata_20220323_NotPublished.Rdata")

