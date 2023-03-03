####################################################################
####################################################################
# --- // Passiflora biflora fertilizer experiment analysis  // --- #
####################################################################
####################################################################


# Colin Morrison and Lauren Hart 
# The University of Texas at Austin 
# Department of Integrative Biology 


#### /////////////////////// ####
#### //// Project GOALS //// ####
#### /////////////////////// ####
# Characterize Passiflora biflora growth & defense traits under different nutrient conditions.
# Quantify whether P. biflora growth and/or defensive abilities are nutrient limited.
# Determine if nutrient ability drives variation if growth or defense of P. biflora.
# Conclude whether P. biflora makes resistance and tolerance trade-offs different nutrient conditions. 
#   This is a test of the resistance-tolerance hypothesis (van der Meijden et al., 1988)

#### /////////////////// ####
#### //// QUESTIONS //// ####
#### /////////////////// ####
# 1.) Are P. biflora growth traits nutrient limited?
#     - HYPOTHESIS: growth will be limited by N and P
#     - ANALYSIS: LME followed by Tukey post-hoc tests
# 2.) Is P. biflora defensive chemical content nutrient limited? 
#     - HYPOTHESIS: cyanogenic glycoside concentration will be limited by N.
#     - ANALYSIS: LME followed by Tukey post-hoc test
# 3.) Are total P. biflora metabolomic richenss and diversity affected by nutrient availability? 
#     - HYPOTHESIS: richness and diveristy will vary with nutrient availability.
#     - ANALYSIS: LME followed by Tukey post-hoc test
# 4.) Does nutrient availability drive shifts in total P. biflora metabolomic profiles? 
#     - HYPOTHESIS: dispersion of groups from different treatments will be distinct in odrination space.
#     - ANALYSIS: PERMANOVA on Bray-Curtis chemical dissimilarity matrices

################################################################################################
################################################################################################
################################################################################################
################################################################################################


getwd()
setwd('~/Desktop/biflora.nutrient.addition')

library(tidyverse)
library(tidyr)
library(readxl)
library(glmm)
library(lme4)
#browseURL('https://www.r-bloggers.com/2017/12/linear-mixed-effect-models-in-r/')
#?lme4
library(nlme)
#?nlme
library(reshape2)
library(dplyr)
library(tidyr)
library(multcomp)
library(MuMIn)
library(plotrix)
library(vegan)


### import data 
data <- read_xlsx("data/databases/total_Pbiflora_dataset_12-16-21.xlsx")
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


#########################################
### subset just final data collection ###
#########################################
t5 <- data %>% filter(timestamp == 5)
View(t5)



########################################################
### --- Question 1.) Effects on Growth Variables --- ###
########################################################
# Growth Variables: SLA, mass_above, mass_below, root.shoot_ratio, internode_length, leaf_number, shoot_length, leafCN_ratio

# QUESTION: Does soil N and P availability affect P. biflora growth traits (length, biomass, SLA, leaf number)?
# HYPOTHESIS: growth poor at low P -> med N/high P-> high N/med --> best at high N and P
# ANALYSES:  1.) LMEs w/ genotype as random effect & weighted variance structure (if data are heteroscedastic)
#            2.) Tukey's post-hoc test for pairwise tests of significance

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

### -- Classic Linear model, no random effects
GLM<-gls(SLA~treatment,method='ML', 
         weights=varIdent(form=~1|as.factor(nitrate)),data=sla.t5)
summary(GLM)
plot(GLM) # residuals normally distributed

### -- LME example testing whether treatment predicted SLA at the final timestamp,
###     genotype modeled as random effect
lmm <- lme(SLA ~ treatment ,random = ~1|genotype, method = "ML", 
           weights=varIdent(form=~1|as.factor(nitrate)),data=sla.t5)
summary(lmm)
plot(lmm) # residuals show homogeneity of variance
###     genotype & individual as nested random effect
lmm2 <- lme(SLA ~ treatment,random = ~1|individual/genotype, method = "ML",
            weights=varIdent(form=~1|as.factor(nitrate)),data=sla.t5)
summary(lmm2)
plot(lmm2) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm2)) 
### These anova() tests will determine if the models are significantly different with respect to goodness-of-fit
# lower the scores, better the model, p-value tells if models are different from each other.
### lmm1 vs lmm2 comparisons
anova(lmm,lmm2) # neither explaining data better (P > 0.05)
### lmm1 vs GLM 
anova(lmm,GLM) # neither explaining data better (P > 0.05)
### lmm12 vs LM 
anova(lmm2,GLM) # models significantly different with respect to goodness-of-fit. 

### So which model?
# GLM has lowest AIC/BIC, but they are very similar. See plots below
anova(lmm,lmm2,GLM)

### --- compare residuals of best fit lme with  GLM (GLM has no random effects, the null)
par(mfrow = c(1,2))
qqnorm(resid(GLM),
       main = "GLM")
abline(0,1, col = "red", lty = 2)
qqnorm(resid(lmm2),
       main = "lmm2")
abline(0,1, col = "red", lty = 2)

### --- Fit final model with restricted maximum likelihood (REML)
# REML estimation is unbiased but does not allow for comparing models with different fixed structures. 
#     So we'll only use the REML estimation on the optimal model.
finalModel <- update(lmm2, .~., method = "REML")
summary(finalModel)

### Now compare with a GLM (no random effects) fit with REML
GLM2 <- update(GLM, .~., method = "REML")
summary(GLM2)

### --- Plot side by side, beta with respective SEs
# We will now contrast our REML-fitted final LME model against a REML-fitted GLM 
#     and determine the impact of incorporating random intercept and slope into the LME, 
#     with respect to SLA and genotype/individual.
# Final GLM
plot(coef(GLM2), xlab = "Fixed Effects", ylab = expression(beta), axes = F,
     pch = 16, col = "black", ylim = c())
stdErrors <- coef(summary(GLM2))[,2]
segments(x0 = 1:9, x1 = 1:9, y0 = coef(GLM2) - stdErrors, y1 = coef(GLM2) + stdErrors,
         col = "black")
axis(2)
abline(h = 0, col = "grey", lty = 2)
axis(1, at = 1:9,
     labels = c("Intercept", "HNLP", "HNMP","LNHP","LNLP",
                "LNMP","MNHP","MNLP","MNMP"), cex.axis = .7)
# Final LME
points(1:9 + .1, fixef(finalModel), pch = 16, col = "red")
stdErrorsLMM <- coef(summary(finalModel))[,2]
segments(x0 = 1:9 + .1, x1 = 1:9 + .1, y0 = fixef(finalModel) - stdErrorsLMM, y1 = fixef(finalModel) + stdErrorsLMM, col = "red")
# Legend
legend("topright", legend = c("GLM","LMM"), text.col = c("black","red"), bty = "n")
# The figure above depicts the estimated beta  from the different fixed effects, 
#     including the intercept, for the GLM (black) and the final LMM (red). 
#     Error bars represent standard errors. 
#     Overall the results are similar for SLA w/ and w/o random effects
# INTERPRETATION: SLA was much higher for plants with low nitrogen regardless of phosphorous,
#                 Including Geno or Geno x Indiv as random effects didn't produce better model.
# WHICH MODEL?: I Let's go with the LME with random effects because it's more  thorough.

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
# --- Leaf area --- #
# ----------------- #
Adata <- read_xlsx("SLA_21area.xlsx")
At5 <- Adata %>% filter(timestamp == 5)

nrow(At5) # 142
At5$leaf_area<-as.numeric(as.character(At5$leaf_area))
At5<-subset(At5, !is.na(leaf_area))
nrow(At5) # [1] 135

### check for normality
Lm1=lm(leaf_area~phosphate*nitrate,data=At5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(leaf_area~nitrate, data = At5) # LOW heterogeneity of variance 
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


# ----------------------- #
# --- Leaf dry weight --- #
# ----------------------- #
Adata <- read_xlsx("SLA_21area.xlsx")
At5 <- Adata %>% filter(timestamp == 5)

nrow(At5) # 142
At5$dry_weight<-as.numeric(as.character(At5$dry_weight))
At5<-subset(At5, !is.na(dry_weight))
nrow(At5) # [1] 135

### check for normality
Lm1=lm(dry_weight~phosphate*nitrate,data=At5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(dry_weight~nitrate, data = At5) # Large heterogeneity of variance 
boxplot(dry_weight~phosphate, data = At5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
At5$nitrate<-as.factor(At5$nitrate)
posthoc<-lme(dry_weight ~ nitrate, random = ~1|individual/genotype.x,method = "REML",
             data=At5) # weights=varIdent(form=~1|nitrate),
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                     Estimate Std. Error z value Pr(>|z|)    
# low - high == 0    -0.022194   0.001601 -13.866  < 1e-10 ***
# medium - high == 0 -0.010143   0.001574  -6.446 1.92e-10 ***
# medium - low == 0   0.012051   0.001592   7.570  < 1e-10 ***


### - compare effect sizes of leaf_area & dry_weight model to which factor drove SLA differences
m1<-lme(dry_weight ~ nitrate,random = ~1|individual/genotype.x,
       method = "REML",data=At5)
summary(m1)
#                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    0.027236   0.001115  24.436  < 2e-16 ***
# nitratelow    -0.022194   0.001594 -13.920  < 2e-16 ***
# nitratemedium -0.010128   0.001559  -6.495 1.55e-09 ***

m2<-lme(leaf_area ~ nitrate,random = ~1|individual/genotype.x,
        method = "REML",data=At5)
summary(m2)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     9.2792     0.4451  20.848  < 2e-16 ***
# nitratelow     -3.0376     0.6367  -4.771 4.79e-06 ***
# nitratemedium  -1.7154     0.6227  -2.755   0.0067 ** 

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


# ------------------ #
# --- mass_total --- #
# ------------------ #
total.t5<-t5[!is.na(t5$mass_total),]

### check for normality
Lm1=lm(mass_total~phosphate*nitrate,data=total.t5)
plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
boxplot(mass_total~nitrate, data = total.t5) #Heterogeneity of variance 
boxplot(mass_total~phosphate, data = total.t5) # LOW heterogeneity of variance 

### -- LME with genotype & individual as nested random effect
lmm <- lme(mass_total ~ nitrate*phosphate,random = ~1|individual/genotype, method = "ML",
           weights=varIdent(form=~1|as.factor(nitrate)),data=total.t5)
plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
finalModel <- update(lmm, .~., method = "REML")
summary(finalModel)
anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   132 755.0806  <.0001
# nitrate               2   132 236.1629  <.0001
# phosphate             2   132   0.3123  0.7323
# nitrate:phosphate     4   132   0.8093  0.5213

r.squaredGLMM(finalModel) 
#            R2m       R2c
# 

### Are the nitrogen levels significantly different from each other?
below.t5$nitrate <- as.factor(below.t5$nitrate)
posthoc<-lme(mass_below ~ nitrate, random = ~1|individual/genotype,
             weights=varIdent(form=~1|as.factor(nitrate)),method = "REML",data=below.t5)
summary(glht(posthoc, linfct=mcp(nitrate="Tukey")))
#                    Estimate Std. Error z value Pr(>|z|)    
#


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
#node.t5<-t5[!is.na(t5$internode_length),]

### check for normality
#Lm1=lm(internode_length~phosphate*nitrate,data=node.t5)
#plot(hist(resid(Lm1))) # normal distribution 
### check for homoscedasticity 
#boxplot(internode_length~nitrate, data = node.t5) # V. Low heterogeneity of variance 
#boxplot(internode_length~phosphate, data = node.t5) # V. low heterogeneity of variance 

### -- LME 
#lmm <- lme(internode_length ~ nitrate*phosphate,
#           random = ~1|individual/genotype, method = "ML",data=node.t5)
#plot(lmm) # residuals show homogeneity of variance

### plot distribution of the random effects
#plot(ranef(lmm)) 

### --- RESULTS:
### Fit final model with REML
#finalModel <- update(lmm, .~., method = "REML")
#summary(finalModel)
#anova(finalModel)
#                   numDF denDF  F-value p-value
# (Intercept)           1   127 1029.6773  <.0001
# nitrate               2   127    0.7336  0.4822
# phosphate             2   127    0.2396  0.7873
# nitrate:phosphate     4   127    1.1326  0.3442

#r.squaredGLMM(finalModel) 
#            R2m       R2c
# [1,] 0.04578079 0.9373141


# ---------------------- #
# --- SOIL N content --- #
# ---------------------- #
soil.t5<-t5[!is.na(t5$soilN),]
nrow(soil.t5) # 18

Lm1=lm(soilN~phosphate*nitrate,data=soil.t5)
plot(hist(resid(Lm1))) # normal distribution 
boxplot(soilN~nitrate, data = soil.t5) # Medium heterogeneity of variance 
boxplot(soilN~phosphate, data = soil.t5) # Low heterogeneity of variance 

### -- LME 
lmm <- lme(soilN ~ nitrate,weights=varIdent(form=~1|as.factor(nitrate)),
           random = ~1|individual/genotype, method = "ML",data=soil.t5)
### --- RESULTS:
summary(lmm)
anova(lmm)
#                   numDF denDF  F-value p-value
# (Intercept)     1    15 3.30457  0.0891
# nitrate         2    15 1.65663  0.2238

View(soil.t5)
stats <- data_summary(soil.t5, varname="soilN", 
                      groupnames=c("nitrate"))
stats
#   nitrate        soilN          sd           se
# 1    high 0.0057857143 0.009710206 0.0036701128
# 2     low 0.0002285714 0.000340168 0.0001285714
# 3  medium 0.0125750000 0.024950000 0.0124750000
### Overall average: 0.0002285714+0.0057857143+0.0125750000/3 = 0.006196429
### Overall SEM:(0.0036701128+0.0001285714+0.0124750000)/3 = 0.005424561

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
### see if pH differed among treatments at end of experiment
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


#########################################################################
### --- Question 2.) Effects on Chemistry - Cyanogenic Glucosides --- ###
#########################################################################
# QUESTION: Does soil N and P availability affect a P. biflora defensive trait (leaf cyanide concentration)?
# HYPOTHESIS: cyanide concentration will increase linearly with N
# ANALYSIS 1: LMEs with genotype as random effect


CG.t5<-t5[!is.na(t5$ug_mg),]
nrow(CG.t5) # [1] 118

### check for CG concentraton normality
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



###########################################################################################
### Question 3.)  Effects on Chemistry - Metabolite Expression - Richness and Diversity ###
###########################################################################################

### --- NPclassifier chemotaxonomic data
npclass=read.table("data/metabolomics/NutrientExp_SiriusQemistreeNPClass_output/ProteoSAFe-NPCLASSIFIER-2ff7c4d9-view_results/NPCLASSIFIER-2ff7c4d9-view_results-main.tsv",
                   header = T, sep = "	", comment.char = "", quote = "")
npclass

### add two columns to NPclass for compounds that have N or P ("containsN" and "containsP" = True) 
npclass = npclass %>% mutate(containsN =
                               case_when(grepl("N",smiles) ~ "True"),
                             containsP =
                               case_when(grepl("P",smiles) ~ "True"))
# fill in cells that for columns that DO NOT have N or P ("containsN" and "containsP" = False) 
npclass = npclass %>% replace_na(list(containsN = 'False',
                                      containsP = 'False'))

### add a column to npclass indicating if a compound is a primary metabolite 
npclass = npclass %>% mutate(metab_type =
                               case_when(pathway_results == "Fatty acids" ~ "primary", 
                                         pathway_results == "Carbohydrates" ~ "primary"))
# fill in the rest of the cells (that have NAs) with secondary metabolite by process of elimination 
npclass = npclass %>% replace_na(list(metab_type = 'secondary')) 

### add a column to npclass indicating if a compound is a cyanogenic glycoside
npclass = npclass %>% mutate(cyanogenic =
                               case_when(grepl("Cyanogenic glycosides",class_results) ~ "True"))
# fill in rest of cells that aren't CGs with NA
npclass = npclass %>% replace_na(list(cyanogenic = 'False'))


### --- Sirius table generated by Qemistree with chemotaxonomic classifications of the metabolites
siri=read.table("data/metabolomics/NutrientExp_SiriusQemistreeNPClass_output/G039_samps_20k_classified_feature_data.tsv",
                header = T, sep = "	", comment.char = "", quote = "")
names(siri)[1] = "id"
siri


### --- combine NPClassifier output with the the "classified_feature_data" Sirius table 
# for-loop to add seven columns from the NP Classifier table to the Sirius table 
for(i in 1:nrow(siri)){
  smiles_i = siri$smiles[i]
  if(smiles_i %in% npclass$smiles){
    siri$class_results[i] = npclass$class_results[which(npclass$smiles == smiles_i)]
    siri$superclass_results[i] = npclass$superclass_results[which(npclass$smiles == smiles_i)]
    siri$pathway_results[i] = npclass$pathway_results[which(npclass$smiles == smiles_i)]
    siri$isglycoside[i] = npclass$isglycoside[which(npclass$smiles == smiles_i)]
    siri$metab_type[i] = npclass$metab_type[which(npclass$smiles == smiles_i)]
    siri$containsN[i] = npclass$containsN[which(npclass$smiles == smiles_i)]
    siri$containsP[i] = npclass$containsP[which(npclass$smiles == smiles_i)]
    siri$cyanogenic[i] = npclass$cyanogenic[which(npclass$smiles == smiles_i)]
  }
}
View(siri)


### --- Feature based molecular network with ion counts (relative abundances) of metabolites in each sample
feat=read.csv("data/metabolomics/NutrientExp_FBMN_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-2014c4a5-download_qza_table_data/quantification_table/quantification_table-00000.csv")
feat
# remove unnecessary characters from sample sames
names(feat) = gsub(pattern = ".mzXML.Peak.area", replacement = "", x = names(feat))  
names(feat) = gsub(pattern = "X", replacement = "", x = names(feat))  
names(feat)
# remove final column with nothing in it
feat=feat[,-154]
dim(feat)
# [1] 2083  153

# fill the matrix of compounds with their ion intensities 
for(i in 1:nrow(siri)){
  labelid = as.character(siri$id[i])
  # featid = as.character(siri$X.featureID[i])
  featid = unlist(strsplit(as.character(siri$X.featureID[i]), split = ","))
  if(length(featid) == 1){
    table = siri$table_number[i]
    if(table == 1){fbmntable = feat}
    # if(table == 1){fbmntable = feat.kanupa44}
    # if(table == 2){fbmntable = feat.titiri42}
    # if(table == 3){fbmntable = feat.tintay25}
    for(j in 4:(ncol(fbmntable)-1)){
      samp = names(fbmntable)[j]
      heat[i,which(names(heat) == samp)] = fbmntable[which(fbmntable$row.ID == featid),j]
    }
    
  }
  if(length(featid) > 1){
    table = unlist(strsplit(as.character(siri$table_number[i]), split = ","))
    for(n in 1:length(featid)){
      # featidn = featid[n]
      if(table[n] == 1){fbmntable = feat}
      # if(table[n] == 1){fbmntable = feat.kanupa44}
      # if(table[n] == 2){fbmntable = feat.titiri42}
      # if(table[n] == 3){fbmntable = feat.tintay25}
      for(j in 4:(ncol(fbmntable)-1)){
        samp = names(fbmntable)[j]
        heat[i,which(names(heat) == samp)] = heat[i,which(names(heat) == samp)] + fbmntable[which(fbmntable$row.ID == featid[n]),j]
      }					
    }
  }	
}
dim(feat)
# 973 153

# remove compounds detected in blanks (3 separate extraction runs w/ 3 unique blanks)
grep("blank", names(heat))# [1]   4   5    117
heat.real = heat[which(heat[,grep("blank1", names(heat))] == 0),]
dim(heat.real)# [1] 794 153
heat.real2 = heat.real[,-c(1:5,117)] # now I can remove the  columns with blank ions and metadata without any information

# take out compounds from siri dataframe that were found in the blanks using heat.real as template
siri.heat = siri[which(siri$id %in% row.names(heat.real2)),] 
dim(siri.heat)# [1] 794  22

# rename X.featureID column row.ID so that siri.heat and feat.real can be merged
names(siri.heat)[names(siri.heat) == 'X.featureID'] <- 'row.ID'
siri.heat=siri.heat[,-c(4:7)] # remove metadata from siri database that contain no information
dim(siri.heat)# [1] 794  19

# now .. take out compounds from feat dataframe that were found in the blanks using siri.heat
feat.real = feat[which(feat$row.ID %in% siri.heat$row.ID),]
dim(feat.real) # [1] 794 153
View(feat.real) # now we have the m/z and retention time data for the final data table
feat.real=feat.real[,-c(4,5,117)] # remove the blanks again (I know this is redundant but it works)


#### --- Merge siri.heat and feat.real so that the m/z and retention times are associated with the ion counts of each sample
master_metab<-merge(feat.real,siri.heat,by="row.ID")
View(master_metab)
#write.csv(master_metab, file="data/metabolomics/master_metab.csv")

# subset the new master metabolimics database into two separate DFs for only samples and only treatment pools
pools_mast_metab=master_metab[,c(1:3,51,129,130,136,139,141,143,145,146,151:168)] # just treatment (nutrient) pools
samps_mast_metab=master_metab[,-c(51,129,130,136,139,141,143,145,146)] # just individual samples (P. biflora plants)

### --- Now melt the dataframe into a long version that shows ion intensities (= presence/absence)
###     of the all compounds in the metabolome by functional chemotypes & biochemical pathway of origin
names(samps_mast_metab)
metab_long <- melt(samps_mast_metab[,c(1,4:141,154:159)], 
                id = c("row.ID","pathway_results","isglycoside","metab_type","containsN","containsP","cyanogenic"))
names(metab_long)=c("row.ID","pathway_results","isglycoside","metab_type","containsN","containsP",
                    "cyanogenic","sample","abundance")

# use the metadata sheet for GNPS  to assign treatments to the samples in a new column of samps_mast_metab
meta=read.table("data/metabolomics/G039_samps_metadata.txt",
                header = T, sep = "	", comment.char = "", quote = "")
head(meta)
meta2=meta[-c(139:141),c(3:7)]
names(meta2)=c("sample","treatment","Nlevel","Plevel","genotype")


### --- Final step in data wrangling is merging metadata and metab_long for statistical analysis
metab_final<-merge(metab_long,meta2,by="sample")
View(metab_final)

### save all these metabolomic data objects for analyses!
save(feat.real,siri.heat,heat.real2,meta,
     master_metab,pools_mast_metab,samps_mast_metab,
     metab_long, metab_final,
     file="data/metabolomics/nutrient_metab.Rdata")



# --------------------------------------------------------- #
# --- Total metabolomic richness per nutrient treatment --- #
# --------------------------------------------------------- #
dats=load(file="data/metabolomics/nutrient_metab.Rdata")
dats 

# Make a table with counts of total unique metabolites/treatment
met_rich <- metab_final %>%
  count(sample, treatment,Nlevel, Plevel, genotype, abundance >0)
names(met_rich)=c("sample","treatment","Nlevel","Plevel","genotype","present","n")

met_rich=met_rich[which(met_rich$present =='TRUE'),names(met_rich)]
nrow(met_rich) # [1] 138 good

# check distribution
Lm=lm(n~Nlevel,data=met_rich)
plot(hist(resid(Lm))) # normal Dist
# check homogeneity of variance
boxplot(n~Nlevel,data=met_rich) # yes heterogeneity of variance 
# model individual and interactive effects  
richness <- lme(n~Nlevel*Plevel,random = ~1|genotype, 
                weights=varIdent(form=~1|as.factor(Nlevel)),
                method = "REML",data=met_rich)
x<-anova(richness) # nitrogen had significant individual effect
# write.csv(x, file="MS/TABLES/TABLE_metab_rich_summary.csv")

# post hoc test for singificant levels
met_rich$Nlevel <- as.factor(met_rich$Nlevel)
posthoc<-lme(n ~ Nlevel, random = ~1|genotype,
             weights=varIdent(form=~1|as.factor(Nlevel)),
             method = "REML",data=met_rich)
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
### low significantly less than medium and high
#                    Estimate Std. Error z value Pr(>|z|)    
# low - high == 0     -57.541      9.934  -5.792   <1e-04 ***
# medium - high == 0  -11.824      7.463  -1.584    0.249    
# medium - low == 0    45.716      9.930   4.604   <1e-04 ***

# ---------------------------------------------------------------------- #
# --- Primary & secondary metabolite richness per nutrient treatment --- #
# ---------------------------------------------------------------------- #
# Make table with counts of unique metabolites/type/treatment
rich_type <- metab_final %>%
  count(sample, treatment,Nlevel, Plevel, genotype, metab_type, abundance >0)
names(rich_type)=c("sample","treatment","Nlevel","Plevel","genotype","metab_type","present","n")

### --- primary metabolites only 
prim_rich=rich_type[which(rich_type$metab_type =='primary' & rich_type$present =='TRUE'),names(rich_type)]
nrow(prim_rich) # [1] 138
# check distribution
Lm=lm(n~Nlevel,data=prim_rich)
plot(hist(resid(Lm))) # normal distribution
# check homogeneity of variance
boxplot(n~Nlevel,data=prim_rich) # very low heterogeneity of variance 
# model individual and interactive effects  
prim_richness <- lme(n~Nlevel*Plevel,random = ~1|genotype, 
                     method = "REML",data=prim_rich)
x<-anova(prim_richness) # nitrogen had significant individual effect
#write.csv(x, file="MS/TABLES/TABLE_primary_rich_summary.csv")

# post hoc test for singificant levels
prim_rich$Nlevel <- as.factor(prim_rich$Nlevel)
posthoc<-lme(n ~ Nlevel, random = ~1|genotype,
             method = "REML",data=prim_rich)
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
### low N significantly less high, medium and high equivalent
#                    Estimate Std. Error z value Pr(>|z|)    
# low - high == 0      -9.143      3.406  -2.685   0.0199 *
# medium - high == 0   -4.259      3.388  -1.257   0.4196  
# medium - low == 0     4.884      3.369   1.450   0.3153 

### --- secondary metabolites only 
second_rich=rich_type[which(rich_type$metab_type =='secondary' & rich_type$present =='TRUE'),names(rich_type)]
nrow(second_rich) # [1] 138
# check distribution
Lm=lm(n~Nlevel,data=second_rich)
plot(hist(resid(Lm))) # nice normal distribution
# check homogeneity of variance
boxplot(n~Nlevel,data=second_rich) # Heterogeneity of variance b/w lows N levels and medN/highN
# model individual and interactive effects  
second_richness <- lme(n~Nlevel*Plevel,random = ~1|genotype, 
                       weights=varIdent(form=~1|as.factor(Nlevel)),
                       method = "REML",data=second_rich)
x<-anova(second_richness) # nitrogen had significant individual effect
#write.csv(x, file="MS/TABLES/TABLE_secondary_rich_summary.csv")

# post hoc test for singificant levels
second_rich$Nlevel <- as.factor(second_rich$Nlevel)
posthoc<-lme(n ~ Nlevel, random = ~1|genotype,
             method = "REML",data=second_rich)
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
### low N significantly less than medium and high
#                    Estimate Std. Error z value Pr(>|z|)    
# low - high == 0     -48.398      7.137  -6.782   <1e-04 ***
# medium - high == 0   -7.565      7.099  -1.066    0.535    
# medium - low == 0    40.832      7.059   5.784   <1e-04 ***

stats <- data_summary(second_rich, varname="n", 
                      groupnames=c("Nlevel"))
stats

# --------------------------------------------------------------------------- #
# ---  richness of compounds with glycoside moeity per nutrient treatment --- #
# --------------------------------------------------------------------------- #
glco_type <- metab_final %>%
  count(sample, treatment,Nlevel, Plevel, genotype, isglycoside, abundance >0)
names(glco_type)=c("sample","treatment","Nlevel","Plevel","genotype","isglycoside","present","n")

glyco_rich=glco_type[which(glco_type$isglycoside =='True' & glco_type$present =='TRUE'),names(glco_type)]
nrow(glyco_rich) # [1] 138
# check distribution
Lm=lm(n~Nlevel,data=glyco_rich)
plot(hist(resid(Lm))) # normal distribution
# check homogeneity of variance
boxplot(n~Nlevel,data=glyco_rich) # yes heterogeneity of variance 
# model individual and interactive effects  
glyco_richness <- lme(n~Nlevel*Plevel,random = ~1|genotype,
                      weights=varIdent(form=~1|as.factor(Nlevel)),
                      method = "REML",data=glyco_rich)
x<-anova(glyco_richness) # nitrogen had significant individual effect
#write.csv(x, file="MS/TABLES/TABLE_isglycoside_rich_summary.csv")

# post hoc test for singificant levels
glyco_rich$Nlevel <- as.factor(glyco_rich$Nlevel)
posthoc<-lme(n ~ Nlevel, random = ~1|genotype,
             method = "REML",data=glyco_rich)
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
### low N significantly less than medium and high
#                    Estimate Std. Error z value Pr(>|z|)    
# low - high == 0      -8.973      1.232  -7.282   <1e-04 ***
# medium - high == 0   -1.927      1.226  -1.572    0.258    
# medium - low == 0     7.047      1.219   5.781   <1e-04 ***


# ---------------------------------------------------------------------- #
# --- N- and P-containing metabolite richness per nutrient treatment --- #
# ---------------------------------------------------------------------- #
### --- Nitrogen-containing metabolites only 
nit_type <- metab_final %>%
  count(sample, treatment,Nlevel, Plevel, genotype, containsN, abundance >0)
names(nit_type)=c("sample","treatment","Nlevel","Plevel","genotype","containsN","present","n")

n_rich=nit_type[which(nit_type$containsN =='True' & nit_type$present =='TRUE'),names(nit_type)]
nrow(n_rich) # [1] 138

# check distribution
Lm=lm(n~Nlevel,data=n_rich)
plot(hist(resid(Lm))) # normal distribution
# check homogeneity of variance
boxplot(n~Nlevel,data=n_rich) # some heterogeneity of variance 
# model individual and interactive effects  
nitro_richness <- lme(n~Nlevel*Plevel,random = ~1|genotype,
                      weights=varIdent(form=~1|as.factor(Nlevel)),
                      method = "REML",data=n_rich)
x<-anova(nitro_richness) # nitrogen had significant individual effect
#write.csv(x, file="MS/TABLES/TABLE_Ncontaining_rich_summary.csv")

# post hoc test for singificant levels
n_rich$Nlevel <- as.factor(n_rich$Nlevel)
posthoc<-lme(n ~ Nlevel, random = ~1|genotype,
             method = "REML",data=n_rich)
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
### low N significantly less than medium and high
# low - high == 0     -32.298      4.645  -6.953   <1e-04 ***
# medium - high == 0   -5.306      4.620  -1.148    0.484    
# medium - low == 0    26.992      4.595   5.875   <1e-04 ***

stats <- data_summary(n_rich, varname="n", 
                      groupnames=c("Nlevel"))
stats

### --- Phosphorus-containing metabolites only 
pho_type <- metab_final %>%
  count(sample,treatment,Nlevel, Plevel, genotype, containsP, abundance >0)
names(pho_type)=c("sample","treatment","Nlevel","Plevel","genotype","containsP","present","n")

p_rich=pho_type[which(pho_type$containsP =='True' & pho_type$present =='TRUE'),names(pho_type)]
nrow(p_rich) # [1] 138

# check distribution
Lm=lm(n~Plevel,data=p_rich)
plot(hist(resid(Lm))) # normal distribution
# check homogeneity of variance
boxplot(n~Plevel,data=p_rich) # very low heterogeneity of variance 
# model individual and interactive effects  
phospho_richness <- lme(n~Nlevel*Plevel,random = ~1|genotype,
                        method = "REML",data=p_rich)
x<-anova(phospho_richness) 
# no patterns of interest. P did not limit number of P-containing compounds.
#write.csv(x, file="MS/TABLES/TABLE_Pcontaining_rich_summary.csv")


# ------------------------------------------------------------ #
# --- Cyanogenic glycoside richness per nutrient treatment --- #
# ------------------------------------------------------------ #
CG_type <- metab_final %>%
  count(sample, treatment,Nlevel, Plevel, genotype, cyanogenic, abundance>0)
names(CG_type)=c("sample","treatment","Nlevel","Plevel","genotype","cyanogenic","present","n")
CG_rich=CG_type[which(CG_type$cyanogenic =='True' & CG_type$present == 'TRUE'),names(CG_type)]
nrow(CG_rich) # [1] 131

# check distribution
Lm=lm(n~Nlevel,data=CG_rich)
plot(hist(resid(Lm))) # normal distribution
# check homogeneity of variance
boxplot(n~Nlevel,data=CG_rich) # some heterogeneity of variance 
# model individual and interactive effects  
CG_richness <- lme(n~Nlevel*Plevel,random = ~1|genotype,
                   weights=varIdent(form=~1|as.factor(Nlevel)),
                   method = "REML",data=CG_rich)
x<-anova(CG_richness) # nitrogen had significant individual effect
#write.csv(x, file="MS/TABLES/TABLE_cyanoGlycoside_rich_summary.csv")

# post hoc test for singificant levels
CG_rich$Nlevel <- as.factor(CG_rich$Nlevel)
posthoc<-lme(n ~ Nlevel, random = ~1|genotype,
             method = "REML",data=CG_rich)
### low N significantly less than medium and high
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
# low - high == 0     -1.5762     0.1977  -7.975  < 1e-04 ***
# medium - high == 0  -0.6375     0.1884  -3.384  0.00201 ** 
# medium - low == 0    0.9387     0.1956   4.798  < 1e-04 ***

# --- what are modal values of unique CGs in each N level?
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Low N 
lowN<-CG_rich[which(CG_rich$Nlevel =='low'),names(CG_rich)]
v<-lowN$n
result <- getmode(v) # [1] 2
# Medium N 
medN<-CG_rich[which(CG_rich$Nlevel =='medium'),names(CG_rich)]
v<-medN$n
result <- getmode(v) # 3 
# High N 
highN<-CG_rich[which(CG_rich$Nlevel =='high'),names(CG_rich)]
v<-highN$n
result <- getmode(v) # 4 

# --- Which cyanogenic glycoside classes are more frequent in different N levels? 
# row.ID 365 = aliphatic CG
# row.ID 1322 = aliphatic CG
# row.ID 412 = aliphatic CG
# row.ID 412 = aromatic CG? Possibly SMC 

CG_type2 <- metab_final %>%
  count(sample, row.ID,  Nlevel, cyanogenic, abundance>0)
names(CG_type2)=c("sample", "row.ID","Nlevel","cyanogenic","present","n")
CG_class=CG_type2[which(CG_type2$cyanogenic =='True' & CG_type2$present == 'TRUE'),names(CG_type2)]
nrow(CG_class) # 
# frequency of each unique CG by N level
CG_count <- CG_class %>%
  count(row.ID,  Nlevel)
CG_count


# ------------------------------------------------------------ #
# --- biochemical pathway diversity per nutrient treatment --- #
# ------------------------------------------------------------ #
names(metab_final)
# Biochemical pathway diversity analyzed with Shannon's Diversity Index.
#   Samples are the 'communities whose diversities will be compared'(= paths_total$sample)
#   Compounds in a sample are analogous to 'species' in a community (= paths_total$pathway).
#   Ion intensities of a given compound are analogous to 'species abundances' (= paths_total$n).  

paths_all <- metab_final %>%
  count(sample, Nlevel, genotype, metab_type, pathway_results, abundance>0)
names(paths_all)=c("sample","Nlevel", "genotype","type", "pathway","present","n")
paths_all=paths_all[which(paths_all$present =='TRUE'),names(paths_all)]

# calculate biochemical pathway diversity for each sample
path_groups=paths_all %>%
  group_by(sample) %>% 
  summarise(val=diversity(n))
path_groups<-as.data.frame(path_groups)
path_groups
# range of H values
range(path_groups$val) # [1] 1.549604 1.920722
# verify function worked by calculating a sample's H by hand - 119
pop = c(8,100,5,15,1,1,61,1,3,18,8)
n = pop
N = sum(pop)
p = n/N
H = -sum(p*log(p))
H # [1] 1.558607 ~ 1.56

# assign treatment level labels back to samples with merge 
path_div<-merge(path_groups,meta2,by="sample")
View(path_div)

### --- all metabolites  
# check distribution
Lm=lm(val~Nlevel,data=path_div)
plot(hist(resid(Lm))) # fairly normal distribution
# check homogeneity of variance
boxplot(val~Nlevel,data=path_div) # some heterogeneity of variance 
# model individual and interactive effects  
all_path_div <- lme(val~Nlevel*Plevel,random = ~1|genotype,
                   weights=varIdent(form=~1|as.factor(Nlevel)),
                   method = "REML",data=path_div)
x<-anova(all_path_div) # nitrogen had significant individual effect
#write.csv(x, file="MS/TABLES/TABLE_biochemPath_total_div_summary.csv")

# post hoc test for singificant levels
path_div$Nlevel <- as.factor(path_div$Nlevel)
posthoc<-lme(val ~ Nlevel, random = ~1|genotype,
             method = "REML",data=path_div)
### low N significantly less than medium and high
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
# low - high == 0    -0.08379    0.01218  -6.879   <1e-04 ***
# medium - high == 0 -0.01778    0.01212  -1.468    0.307    
# medium - low == 0   0.06601    0.01205   5.478   <1e-04 ***


### --- primary metabolites only 
paths_prim=paths_all[which(paths_all$present =='TRUE' & paths_all$type =='primary'),names(paths_all)]
# calculate biochemical pathway diversity for each sample
prim_groups=paths_prim %>%
  group_by(sample) %>% 
  summarise(val=diversity(n))
prim_groups<-as.data.frame(prim_groups)
prim_groups
# range of H values
range(prim_groups$val) # [1] 0.1461447 0.4055999
# assign treatment level labels back to samples with merge 
prim_div<-merge(prim_groups,meta2,by="sample")
View(prim_div)

# check distribution
Lm=lm(val~Nlevel,data=prim_div)
plot(hist(resid(Lm))) # fairly normal distribution
# check homogeneity of variance
boxplot(val~Nlevel,data=prim_div) # some heterogeneity of variance 
# model individual and interactive effects  
prim_path_div <- lme(val~Nlevel*Plevel,random = ~1|genotype,
                    weights=varIdent(form=~1|as.factor(Nlevel)),
                    method = "REML",data=prim_div)
x<-anova(prim_path_div) # no significant effects
#write.csv(x, file="MS/TABLES/TABLE_biochemPath_primary_div_summary.csv")


### --- secondary metabolites only 
paths_second=paths_all[which(paths_all$present =='TRUE' & paths_all$type =='secondary'),names(paths_all)]
# calculate biochemical pathway diversity for each sample
second_groups=paths_second %>%
  group_by(sample) %>% 
  summarise(val=diversity(n))
second_groups<-as.data.frame(second_groups)
second_groups
# range of H values
range(second_groups$val) # [1] 1.349584 1.722246
# assign treatment level labels back to samples with merge 
second_div<-merge(second_groups,meta2,by="sample")
View(second_div)

# check distribution
Lm=lm(val~Nlevel,data=second_div)
plot(hist(resid(Lm))) # fairly normal distribution
# check homogeneity of variance
boxplot(val~Nlevel,data=second_div) # some heterogeneity of variance 
# model individual and interactive effects  
second_path_div <- lme(val~Nlevel*Plevel,random = ~1|genotype,
                     weights=varIdent(form=~1|as.factor(Nlevel)),
                     method = "REML",data=second_div)
x<-anova(second_path_div) # nitrogen had significant individual effect
#write.csv(x, file="MS/TABLES/TABLE_biochemPath_secondary_div_summary.csv")

# post hoc test for singificant levels
second_div$Nlevel <- as.factor(second_div$Nlevel)
posthoc<-lme(val ~ Nlevel, random = ~1|genotype,
             method = "REML",data=second_div)
### low N significantly less than medium and high
summary(glht(posthoc, linfct=mcp(Nlevel="Tukey")))
# low - high == 0    -0.08254    0.01149  -7.186   <0.001 ***
# medium - high == 0 -0.02561    0.01143  -2.242   0.0645 .  
# medium - low == 0   0.05693    0.01136   5.010   <0.001 ***


### --- save all these NPclassifier group objects for plotting
save(metab_final,meta2,
     met_rich,rich_type,prim_rich,second_rich,
     glyco_rich,n_rich,p_rich,CG_rich,CG_count,
     path_div,prim_div,second_div,
     file="data/metabolomics/metabolomics_data.Rdata")



#########################################################################################
### Question 4.)  Effects on Shifts in Metabolomic Expression across Ordination Space ###
#########################################################################################
dats=load(file="CSCS_biflora_nutrient-primary_20221026.RData")
dats

### Remove species pools from matrix
dim(braydist.bin) # [1] 147 147
braydist.bin.samps = braydist.bin[,-grep("_pool",names(braydist.bin))] # remove columns with pools
row.names(braydist.bin.samps)
braydist.bin.samps = braydist.bin.samps[-c(48,126,127,133,136,138,140,142,143),] # remove rows with pools 
dim(braydist.bin.samps) # [1] 138 138
### calculate NMDS axes using vegan
library(vegan)
nmds.samps2 = metaMDS(as.dist(braydist.bin.samps), k=5)
nmds.samps.points2 = as.data.frame(nmds.samps2$points)
### add column with sample names to NMDS points df for merging with meta data file
nmds.samps.points2$ATTRIBUTE_SampleNumber=row.names(nmds.samps.points2)
### merge with treatment categories from meta data file
nmds.ready2<-merge(nmds.samps.points2,meta[,c(3:6)],by="ATTRIBUTE_SampleNumber")

### --- PERMANOVA with Adonis
lowVmed<-nmds.ready2[-which(nmds.ready2$ATTRIBUTE_Nitrogen == "high"),names(nmds.ready2)]
lowVhigh<-nmds.ready2[-which(nmds.ready2$ATTRIBUTE_Nitrogen == "medium"),names(nmds.ready2)]
medVhigh<-nmds.ready2[-which(nmds.ready2$ATTRIBUTE_Nitrogen == "low"),names(nmds.ready2)]

axes<-data.frame(lowVmed[,2:6])
manova <- adonis(axes ~ ATTRIBUTE_Nitrogen, method = "euclidean",data=lowVmed,
                 permutations = 9999, sqrt.dist = TRUE)
manova

### --- WHOLE plant metabolomes
### Low versus medium N 
#                     Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1    0.2996 0.299645  8.5183 0.0856  1e-04 ***
# Residuals          91    3.2011 0.035176         0.9144           
# Total              92    3.5007                  1.0000                    

### Low versus high N 
#                    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1    0.6160 0.61603  18.249 0.17016  1e-04 ***
# Residuals          89    3.0043 0.03376         0.82984           
# Total              90    3.6204                 1.00000                      

### Medium versus high N
#                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1   0.11246 0.112462  4.4577 0.04719 0.0018 **
# Residuals          90   2.27058 0.025229         0.95281          
# Total              91   2.38304                  1.00000         


### --- PRIMARY metabolites only
### Low versus medium N 
#                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1    0.2330 0.232950  5.9058 0.06094  2e-04 ***
# Residuals          91    3.5894 0.039444         0.93906           
# Total              92    3.8224                  1.00000         

### Low versus high N 
#                    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1    0.5518 0.55183   14.39 0.13918  1e-04 ***
# Residuals          89    3.4130 0.03835         0.86082           
# Total              90    3.9648                 1.00000 

### Medium versus high N
#                    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1    0.1491 0.14909  4.4227 0.04684 0.0024 **
# Residuals          90    3.0339 0.03371         0.95316          
# Total              91    3.1830                 1.00000          


### --- SECONDARY metabolites only
### Low versus medium N 
#                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1    0.3979 0.39786  10.197 0.10076  1e-04 ***
# Residuals          91    3.5507 0.03902         0.89924           
# Total              92    3.9486                 1.00000        

### Low versus high N 
#                    Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1    0.7245 0.72449  19.491 0.17966  1e-04 ***
# Residuals          89    3.3081 0.03717         0.82034           
# Total              90    4.0326                 1.00000                  

### Medium versus high N
#                   Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# ATTRIBUTE_Nitrogen  1   0.10234 0.102341  4.2539 0.04513 0.0018 **
# Residuals          90   2.16522 0.024058         0.95487          
# Total              91   2.26756                  1.00000    



### //////////////////////////////////////////////////////////////////////// ###
### /// TABLE S. NMDS Ordination Figure of Total Metabolomic Similarity /// ###
### //////////////////////////////////////////////////////////////////////// ###
dats=load(file="data/metabolomics/nutrient_metab.Rdata")
dats 

View(pools_mast_metab)
dim(pools_mast_metab) # [1] 794  30

### remove rows with 0 ion count in each sample (means that compound wasn't found in any of the the pools)
pools_mast_metab<-pools_mast_metab[!duplicated(pools_mast_metab[c(4:12)]),]

dim(pools_mast_metab) # [1] 593  30

write.csv(pools_mast_metab, file="MS/TABLES/TABLE_S5_biflora_metabolome.csv")



###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################