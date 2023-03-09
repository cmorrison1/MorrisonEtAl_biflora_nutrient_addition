####################################################################
####################################################################
# --- // Passiflora biflora fertilizer experiment analysis  // --- #
####################################################################
####################################################################

# Colin Morrison and Lauren Hart 
# The University of Texas at Austin 
# Department of Integrative Biology 


setwd('~/Documents/CRM manuscripts AND essays/Biflora.nutrient.addition')

library(tidyverse)
library(readxl)
library(lme4)
library(nlme)
library(reshape2)
library(dplyr)
library(tidyr)
library(multcomp)
library(MuMIn)
library(plotrix)


### import data 
data <- read_xlsx("total_Pbiflora_dataset_12-16-21.xlsx")
data=data[,-1]
View(data)
### add new column for phosphate level 
data<- data %>%
  mutate(phosphate = case_when(
    endsWith(treatment, "LP") ~ "low",
    endsWith(treatment, "MP") ~ "medium",
    endsWith(treatment, "HP") ~ "high"
  ))
### add new column for nitrate level 
data<- data %>%
  mutate(nitrate = case_when(
    startsWith(treatment, "LN") ~ "low",
    startsWith(treatment, "MN") ~ "medium",
    startsWith(treatment, "HN") ~ "high"
  ))

data<-data[,c(1:6,46:47,7:45)] # reorder the columns so all predictors up front
names(data)


### make summary data table with mean, sd, and se for the plot
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- plyr::rename(data_sum, c("mean" = varname))
  return(data_sum)
}



#####################################################
### --- Goal 3b.) Effects on Growth Variables --- ###
#####################################################
# Growth Variables: SLA, mass_above, mass_below, root.shoot_ratio, internode_length, leaf_number, shoot_length, leafCN_ratio

# QUESTION: Does soil N and P availability affect P. biflora growth traits (length, biomass, SLA, leaf number)?
# HYPOTHESIS: growth poor at low P -> med N/high P-> high N/med --> best at high N and P
# ANALYSES:  1.) LMEs w/ genotype as random effect & weighted variance structure (if data are heteroscedastic)
#            2.) Tukey's post-hoc test for pairwise tests of significance



### subset just final data collection 
t5 <- data %>% filter(timestamp == 5)
View(t5)


# ----------- #
# --- SLA --- #
# ----------- #
### NOTE: Higher SLA = thinner leaves; lower SLA = thicker leaves
nrow(t5) # 142
sla.t5<-t5[!is.na(t5$SLA),]
nrow(sla.t5) # [1] 134

### check for normality
Lm1=lm(SLA~phosphate*nitrate,data=sla.t5)
summary(Lm1)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(SLA~nitrate, data = sla.t5) # Large heterogeneity of variance -> weight the predictors
boxplot(SLA~phosphate, data = sla.t5) # Large heterogeneity of variance -> weight the predictors


### -- LME example testing whether treatment predicted SLA at the final timestamp,
###     genotype & individual as nested random effect
lmm2 <- lme(SLA ~ treatment,random = ~1|individual/genotype, method = "ML",
            weights=varIdent(form=~1|as.factor(nitrate)),data=sla.t5)
summary(lmm2)
plot(lmm2) # residuals show homogeneity of variance

### --- Fit final model with restricted maximum likelihood (REML)
# REML estimation is unbiased but does not allow for comparing models with different fixed structures. 
#     So we'll only use the REML estimation on the optimal model.
finalModel <- update(lmm2, .~., method = "REML")
summary(finalModel)

### --- Results:
finalModel2 <- lme(SLA ~ nitrate*phosphate,random = ~1|individual/genotype, method = "REML",
                   weights=varIdent(form=~1|as.factor(nitrate)),data=sla.t5)
summary(finalModel2)
anova(finalModel2)
#                   numDF denDF  F-value p-value
# (Intercept)           1   125 356.9706  <.0001
# nitrate               2   125  43.4458  <.0001
# phosphate             2   125   0.1606  0.8518
# nitrate:phosphate     4   125   0.5201  0.7211

r.squaredGLMM(finalModel2) 
#            R2m       R2c
# [1,] 0.7600643 0.9656794   ### WOW!

### make data table with mean, sd, and se for summary statistics
stats <- data_summary(t5, varname="SLA", 
                      groupnames=c("treatment"))
stats


# ----------------- #
# --- Leaf Area --- #
# ----------------- #
Adata <- read_xlsx("SLA_21area.xlsx")
At5 <- Adata %>% filter(timestamp == 5)

nrow(At5) # 142
At5$leaf_area<-as.numeric(as.character(At5$leaf_area))
At5<-subset(At5, !is.na(leaf_area))
nrow(At5) # [1] 142

### check for normality
Lm1=lm(leaf_area~phosphate*nitrate,data=At5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(leaf_area~nitrate, data = At5) # Some heterogeneity of variance 
boxplot(leaf_area~phosphate, data = At5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(leaf_area ~ nitrate*phosphate,random = ~1|individual/genotype.x, method = "ML",
           weights=varIdent(form=~1|as.factor(nitrate)),data=At5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# 
r.squaredGLMM(finalModel2) 
#            R2m       R2c
# (Intercept)           1   125 1063.0415  <.0001
# nitrate               2   125   10.6081  0.0001
# phosphate             2   125    1.6258  0.2009
# nitrate:phosphate     4   125    1.1141  0.3529


### Are the nitrogen levels significantly different from each other?
At5$nitrate <- as.factor(At5$nitrate)
posthoc<-lme(leaf_area ~ nitrate, random = ~1|individual/genotype.x,
             method = "REML",data=At5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|)
# low - high == 0     -3.0376     0.6367  -4.771   <0.001 ***
# medium - high == 0  -1.7154     0.6227  -2.755   0.0162 *  
# medium - low == 0    1.3222     0.6301   2.098   0.0901 .  



# ------------------ #
# --- mass_above --- #
# ------------------ #
nrow(t5) # 142
above.t5<-t5[!is.na(t5$mass_above),]
nrow(above.t5) # [1] 141

### check for normality
Lm1=lm(mass_above~phosphate*nitrate,data=above.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(mass_above~nitrate, data = above.t5) # Large heterogeneity of variance 
boxplot(mass_above~phosphate, data = above.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(mass_above ~ nitrate*phosphate,random = ~1|individual/genotype, method = "ML",
           weights=varIdent(form=~1|as.factor(nitrate)),data=above.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   132 681.6140  <.0001
# nitrate               2   132 267.4674  <.0001
# phosphate             2   132   0.0604  0.9414
# nitrate:phosphate     4   132   0.9322  0.4474
r.squaredGLMM(finalModel2) 
#            R2m       R2c
# [1,] 0.5320848 0.5321087

### Are the nitrogen levels significantly different from each other?
above.t5$nitrate <- as.factor(above.t5$nitrate)
posthoc<-lme(mass_above ~ nitrate, random = ~1|individual/genotype,
             method = "REML",data=above.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|)    
# low - high == 0    -1.59431    0.07927 -20.112   <1e-10 ***
# medium - high == 0 -0.94516    0.08007 -11.804   <1e-10 ***
# medium - low == 0   0.64915    0.07839   8.281   <1e-10 ***

# ------------------ #
# --- mass_below --- #
# ------------------ #
below.t5<-t5[!is.na(t5$mass_below),]

### check for normality
Lm1=lm(mass_below~phosphate*nitrate,data=below.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(mass_below~nitrate, data = below.t5) # Some heterogeneity of variance 
boxplot(mass_below~phosphate, data = below.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(mass_below ~ nitrate*phosphate,random = ~1|individual/genotype, method = "ML",
           weights=varIdent(form=~1|as.factor(nitrate)),data=below.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   132 485.3544  <.0001
# nitrate               2   132  99.8538  <.0001
# phosphate             2   132   1.0428  0.3553
# nitrate:phosphate     4   132   0.8206  0.5143
r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.2655169 0.2704826

### Are the nitrogen levels significantly different from each other?
below.t5$nitrate <- as.factor(below.t5$nitrate)
posthoc<-lme(mass_below ~ nitrate, random = ~1|individual/genotype,
             weights=varIdent(form=~1|as.factor(nitrate)),method = "REML",data=below.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|)    
#low - high == 0    -0.45969    0.04809  -9.558  < 1e-10 ***
# medium - high == 0 -0.29666    0.04981  -5.955 6.76e-09 ***
# medium - low == 0   0.16303    0.01492  10.925  < 1e-10 ***



# ------------------------ #
# --- root.shoot_ratio --- #
# ------------------------ #
rs.t5<-t5[!is.na(t5$root.shoot_ratio),]

### check for normality
Lm1=lm(root.shoot_ratio~phosphate*nitrate,data=rs.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(root.shoot_ratio~nitrate, data = rs.t5) # Some heterogeneity of variance 
boxplot(root.shoot_ratio~phosphate, data = rs.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(root.shoot_ratio ~ nitrate*phosphate,random = ~1|individual/genotype, method = "ML",
           weights=varIdent(form=~1|as.factor(nitrate)),data=rs.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   132 1208.9514  <.0001
# nitrate               2   132   10.3175  0.0001
# phosphate             2   132    0.0398  0.9610
# nitrate:phosphate     4   132    1.6245  0.1718
r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.3137698 0.6081427

### Are the nitrogen levels significantly different from each other?
rs.t5$nitrate <- as.factor(rs.t5$nitrate)
posthoc<-lme(root.shoot_ratio ~ nitrate, random = ~1|individual/genotype,
             method = "REML",data=rs.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|)    
# low - high == 0     0.13507    0.03041   4.441 3.11e-05 ***
# medium - high == 0 -0.01138    0.03072  -0.370    0.927    
# medium - low == 0  -0.14644    0.03007  -4.870  < 1e-05 ***


# ------------------- #
# --- leaf_number --- #
# ------------------- #
leaf.t5<-t5[!is.na(t5$leaf_number),]

### check for normality
Lm1=lm(leaf_number~phosphate*nitrate,data=leaf.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(leaf_number~nitrate, data = leaf.t5) # BIG heterogeneity of variance 
boxplot(leaf_number~phosphate, data = leaf.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(leaf_number ~ nitrate*phosphate,random = ~1|individual/genotype, method = "ML",
           weights=varIdent(form=~1|as.factor(nitrate)),data=leaf.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   130 1041.1643  <.0001
# nitrate               2   130  275.9857  <.0001
# phosphate             2   130    0.9873  0.3754
# nitrate:phosphate     4   130    0.4239  0.7911

r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.6406281 0.6467881

### Are the nitrogen levels significantly different from each other?
leaf.t5$nitrate <- as.factor(leaf.t5$nitrate)
posthoc<-lme(leaf_number ~ nitrate, random = ~1|individual/genotype,
             method = "REML",data=leaf.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|) 
# low - high == 0     -22.978      1.118 -20.549   <2e-16 ***
# medium - high == 0   -9.539      1.135  -8.402   <2e-16 ***
# medium - low == 0    13.440      1.105  12.159   <2e-16 ***

# -------------------- #
# --- shoot_length --- #
# -------------------- #
shoot.t5<-t5[!is.na(t5$shoot_length),]

### check for normality
Lm1=lm(shoot_length~phosphate*nitrate,data=shoot.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(shoot_length~nitrate, data = shoot.t5) # BIG heterogeneity of variance 
boxplot(shoot_length~phosphate, data = shoot.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(shoot_length ~ nitrate*phosphate,random = ~1|individual/genotype, method = "ML",
           weights=varIdent(form=~1|as.factor(nitrate)),data=shoot.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   130 1077.0613  <.0001
# nitrate               2   130  207.4413  <.0001
# phosphate             2   130    0.4282  0.6526
# nitrate:phosphate     4   130    1.7429  0.1444

r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.5050392 0.51616

### Are the nitrogen levels significantly different from each other?
shoot.t5$nitrate <- as.factor(shoot.t5$nitrate)
posthoc<-lme(shoot_length ~ nitrate, random = ~1|individual/genotype,
             method = "REML",data=shoot.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|) 
# low - high == 0    -153.303      8.600 -17.825  < 2e-16 ***
# medium - high == 0  -72.199      8.732  -8.268 3.33e-16 ***
# medium - low == 0    81.104      8.501   9.540  < 2e-16 ***

### make data table with mean, sd, and se for summary statistics
stats <- data_summary(shoot.t5, varname="shoot_length", 
                      groupnames=c("treatment"))
stats


# ------------------------ #
# --- Internode length --- #
# ------------------------ #
node.t5<-t5[!is.na(t5$internode_length),]

### check for normality
Lm1=lm(internode_length~phosphate*nitrate,data=node.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(internode_length~nitrate, data = node.t5) # V. Low heterogeneity of variance 
boxplot(internode_length~phosphate, data = node.t5) # V. low heterogeneity of variance 

### -- LME 
lmm <- lme(internode_length ~ nitrate*phosphate,
           random = ~1|individual/genotype, method = "ML",data=node.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   127 1029.6773  <.0001
# nitrate               2   127    0.7336  0.4822
# phosphate             2   127    0.2396  0.7873
# nitrate:phosphate     4   127    1.1326  0.3442

r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.04578079 0.9373141


# ---------------------- #
# --- leaf N content --- #
# ---------------------- #
N.t5<-t5[!is.na(t5$leafN),]
nrow(N.t5) # [1] 106

### check for normality
Lm1=lm(leafN~phosphate*nitrate,data=N.t5) # one huge outlier
plot(hist(resid(Lm1))) # normal distribution 
N.t5=N.t5[-86,] # remove the outlier to conform to normality 
### check for homoscedasticity 
boxplot(leafN~nitrate, data = N.t5) # Medium heterogeneity of variance 
boxplot(leafN~phosphate, data = N.t5) # V. Low heterogeneity of variance 

### -- LME 
lmm <- lme(leafN ~ nitrate*phosphate,weights=varIdent(form=~1|as.factor(nitrate)),
           random = ~1|individual/genotype, method = "ML",data=N.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1    96 1618.7571  <.0001
# nitrate               2    96   20.1035  <.0001
# phosphate             2    96    0.3544  0.7025
# nitrate:phosphate     4    96    0.7353  0.5701

r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.2672836 0.7407882

N.t5$nitrate <- as.factor(N.t5$nitrate)
posthoc<-lme(leafN ~ nitrate, random = ~1|individual/genotype,
             weights=varIdent(form=~1|as.factor(nitrate)), method = "REML",data=N.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|) 
# low - high == 0    -0.84556    0.20653  -4.094 0.000117 ***
# medium - high == 0 -0.90523    0.14602  -6.199  < 1e-04 ***
# medium - low == 0  -0.05967    0.19644  -0.304 0.949660  


# ---------------------- #
# --- leaf C content --- #
# ---------------------- #
C.t5<-t5[!is.na(t5$leafC),]
C.t5=C.t5[-86,] # N.t5=N.t5[-86,] # remove the outlier to conform to normality (ESI processing error)

### check for normality
Lm1=lm(leafC~phosphate*nitrate,data=C.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(leafC~nitrate, data = C.t5) # Low heterogeneity of variance 
boxplot(leafC~phosphate, data = C.t5) # Low heterogeneity of variance 

### -- LME 
lmm <- lme(leafC ~ nitrate*phosphate, random = ~1|individual/genotype, 
           method = "ML",data=C.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1    95 5993.586  <.0001
# nitrate               2    95    0.762  0.4695
# phosphate             2    95    1.096  0.3383
# nitrate:phosphate     4    95    0.870  0.4848



# -------------------- #
# --- leafCN_ratio --- #
# -------------------- #
CN.t5<-t5[!is.na(t5$leafCN_ratio),]
nrow(CN.t5) # [1] 106

### check for normality
Lm1=lm(leafCN_ratio~phosphate*nitrate,data=CN.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(leafCN_ratio~nitrate, data = CN.t5) # Some heterogeneity of variance 
boxplot(leafCN_ratio~phosphate, data = CN.t5) # V. low heterogeneity of variance 

### -- LME 
lmm <- lme(leafCN_ratio ~ nitrate*phosphate,weights=varIdent(form=~1|as.factor(nitrate)),
           random = ~1|individual/genotype, method = "ML",data=CN.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1    97 1874.7712  <.0001
# nitrate               2    97   21.2158  <.0001
# phosphate             2    97    1.0595  0.3506
# nitrate:phosphate     4    97    0.8187  0.5163

r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.4200063 0.9220178

CN.t5$nitrate <- as.factor(CN.t5$nitrate)
posthoc<-lme(leafCN_ratio ~ nitrate, random = ~1|individual/genotype,
             weights=varIdent(form=~1|as.factor(nitrate)), method = "REML",data=CN.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|) 
# low - high == 0      4.1187     1.5975   2.578   0.0248 *  
# medium - high == 0   5.1955     0.8179   6.352   <0.001 ***
# medium - low == 0    1.0768     1.6653   0.647   0.7863  

### make data table with mean, sd, and se for summary statistics
stats <- data_summary(CN.t5, varname="leafCN_ratio", 
                      groupnames=c("treatment"))
stats

# --------------- #
# --- soil pH --- #
# --------------- #
nrow(t5) # 142
pH.t5<-t5[!is.na(t5$pH),]
nrow(pH.t5) # [1] 140

### check for normality
Lm1=lm(pH~phosphate*nitrate,data=pH.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(pH~nitrate, data = pH.t5) # LOW heterogeneity of variance 
boxplot(pH~phosphate, data = pH.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(pH ~ treatment,random = ~1|individual/genotype, method = "ML",data=pH.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   131 108937.27  <.0001
# nitrate               2   131     0.84  0.4329
# phosphate             2   131     1.66  0.1947
# nitrate:phosphate     4   131     0.69  0.6024
r.squaredGLMM(finalModel2) 
#            R2m       R2c
# [1,] 0.05278263 0.937774


########################################################
### --- Goal 3c.) Effects on Defensive Chemistry --- ###
########################################################
# QUESTION: Does soil N and P availability affect a P. biflora defensive trait (leaf cyanide concentration)?
# HYPOTHESIS: cyanide concentration will increase linearly with N
# ANALYSIS: LMEs with genotype as random effect


CG.t5<-t5[!is.na(t5$ug_mg),]
nrow(CG.t5) # [1] 118

### check for normality
Lm1=lm(ug_mg~phosphate*nitrate,data=CG.t5)
plot(hist(resid(Lm1))) # huge outlier on right
CG.t5=CG.t5[-93,] # remove outlier to conform to norm assumption (it doesn't change the LME or posthoc results)
### check for homoscedasticity 
par(mfrow=c(1,2))
boxplot(ug_mg~nitrate, data = CG.t5) # Some heterogeneity of variance 
boxplot(ug_mg~phosphate, data = CG.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(ug_mg ~ nitrate*phosphate,random = ~1|individual/genotype, 
           weights=varIdent(form=~1|as.factor(nitrate)),method = "ML",data=CG.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   108 37.14183  <.0001
# nitrate               2   108 57.80644  <.0001
# phosphate             2   108  0.33841  0.7137
# nitrate:phosphate     4   108  0.56854  0.6860
r.squaredGLMM(finalModel2) 
#            R2m       R2c
# [1,] 0.5320848 0.5321087

CG.t5$nitrate <- as.factor(CG.t5$nitrate)
posthoc<-lme(ug_mg ~ nitrate, random = ~1|individual/genotype,
             weights=varIdent(form=~1|as.factor(nitrate)), method = "REML",data=CG.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|) 
#  low - high == 0    -0.56921    0.07743  -7.352   <1e-04 ***
# medium - high == 0 -0.21026    0.08714  -2.413   0.0383 *  
# medium - low == 0   0.35895    0.04316   8.317   <1e-04 ***

### make data table with mean, sd, and se for summary statistics
stats <- data_summary(CG.t5, varname="ug_mg", 
                      groupnames=c("treatment"))
stats


####################################################################
####################################################################
####################################################################
####################################################################