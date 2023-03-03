################################################################
################################################################
# --- // Passiflora biflora fertilizer experiment PLOTS // --- #
################################################################
################################################################



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
#       This is a test of the resistance-tolerance hypothesis (van der Meijden et al., 1988)


library(tidyverse)
library(readxl)
library(lme4)
library(nlme)
library(reshape2)
library(dplyr)
library(tidyr)
library(MuMIn)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(vegan)

getwd()
setwd('~/Desktop/biflora.nutrient.addition')

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


##############################################################
##############################################################
##############################################################
### --- Goal 1.) Nitrogen effect on chemistry in Field --- ###
##############################################################
##############################################################
##############################################################
# QUESTION: Does soil N and P availability affect concentration of a N-based leaf secondary defense?
# HYPOTHESIS: There will be a correlation b/w soil NH4-NO3 and P. biflora leaf HCN 
# PLOTS: scatter plots with 95% confidence intervals 

data=load(file="data/hcn.soil.nitrogen_LS/LS_field_N-HCN.Rdata")
# [1] "fielDat"


### --- Plot it
LS<-ggscatter(fielDat, x = "totalN.log", y = "hcn.log2",
              add = "reg.line",                     # Add regression line
              conf.int = TRUE,                      # Add confidence interval
              size = 3,
              add.params = list(color = "blue",
                                fill = "lightgray")) +
  annotate(geom="text", x=1.63, y=0.18, label=expression(paste(italic('P'), " = 0.02" )),
           size=6) + 
  annotate(geom="text", x=1.63, y=0.17, label=expression(paste('R'^2, ' = 0.13')),
           size=6) + 
  xlab(expression(paste('Log'[10], ' soil ammonium nitrate (ug/mg)'))) +
  ylab(expression(paste('Log'[10], ' leaf hydrogen cyanide (ug/mg)'))) +
  theme(panel.grid.major.x = element_line(colour = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=-1),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )
LS
#

dev.new(
  title = "LS HCN-soil N",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
LS



####################################################
####################################################
####################################################
### --- Goal 2.) Effects on Growth Variables --- ###
####################################################
####################################################
####################################################
# QUESTION: Does soil N and P availability affect P. biflora growth traits (length, biomass, SLA, leaf number)?
# HYPOTHESIS: growth poor at low P -> med N/high P-> high N/med --> best at high N and P
# PLOTS: errorplots - average values plotted with standard error bars on both sides

# import data 
data <- read_xlsx("data/databases/total_Pbiflora_dataset_12-16-21.xlsx")
data=data[,-1]
View(data)

### --- add new column for phosphate level 
#data<- data %>%
#  mutate(phosphate = case_when(
#    endsWith(treatment, "LP") ~ "low",
#    endsWith(treatment, "MP") ~ "medium",
#    endsWith(treatment, "HP") ~ "high"
#  ))
###  --- add new column for nitrate level 
#data<- data %>%
#  mutate(nitrate = case_when(
#    startsWith(treatment, "LN") ~ "low",
#    startsWith(treatment, "MN") ~ "medium",
#    startsWith(treatment, "HN") ~ "high"
#  ))



### subset just final data collection 
t5 <- data %>% filter(timestamp == 5)
# View(t5)


### /////////////////////////////////////////////////////////////////// ###
### /// FIGURE 2. Effects on P.biflora Growth and Functional Traits /// ###
### /////////////////////////////////////////////////////////////////// ###

# --------------------------- #
# --- FIG. 2A: mass_total --- #
# --------------------------- #
total.t5<-t5[!is.na(t5$mass_total),]
total.t5$treatment <- factor(total.t5$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

total<-ggerrorplot(total.t5, x = "treatment", y = "mass_total", 
                   desc_stat = "mean_se", color = "black",
                   size=0.5,width = 2.5) +
  ylab('Total biomass (g)') +
  xlab('Nutrient treatment') + 
  scale_y_continuous(breaks =seq(0,3.0,0.5)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "total biomass",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
total


# --------------------------------- #
# --- FIG. 2B: root.shoot_ratio --- #
# --------------------------------- #
rs.t5<-t5[!is.na(t5$root.shoot_ratio),]
rs.t5$treatment <- factor(rs.t5$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

RS<-ggerrorplot(rs.t5, x = "treatment", y = "root.shoot_ratio", 
                   desc_stat = "mean_se", color = "black",
                   size=0.5,width = 2.5) +
  ylab('Root/shoot ratio') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "root:shoot",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
RS

# ---------------------------- #
# --- FIG. 2C: leaf_number --- #
# ---------------------------- #
leaf.t5<-t5[!is.na(t5$leaf_number),]
leaf.t5$treatment <- factor(leaf.t5$treatment,
                          levels = c("LNLP", "LNMP", "LNHP",
                                     "MNLP", "MNMP", "MNHP",
                                     "HNLP", "HNMP", "HNHP"))

leaf<-ggerrorplot(leaf.t5, x = "treatment", y = "leaf_number", 
                desc_stat = "mean_se", color = "black",
                size=0.5,width = 2.5) +
  ylab('Leaf number') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "leaves",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
leaf

# ----------------------------- #
# --- FIG. 2D: shoot_length --- #
# ----------------------------- #
shoot.t5<-t5[!is.na(t5$shoot_length),]
shoot.t5$treatment <- factor(shoot.t5$treatment,
                            levels = c("LNLP", "LNMP", "LNHP",
                                       "MNLP", "MNMP", "MNHP",
                                       "HNLP", "HNMP", "HNHP"))

shoot<-ggerrorplot(shoot.t5, x = "treatment", y = "shoot_length", 
                  desc_stat = "mean_se", color = "black",
                  size=0.5,width = 2.5) +
  ylab('Shoot length (cm)') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "shoots",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
shoot


# -------------------- #
# --- FIG. 2E: SLA --- #
# -------------------- #
sla.t5<-t5[!is.na(t5$SLA),]
sla.t5$treatment <- factor(sla.t5$treatment,
                           levels = c("LNLP", "LNMP", "LNHP",
                                      "MNLP", "MNMP", "MNHP",
                                      "HNLP", "HNMP", "HNHP"))

SLA<-ggerrorplot(sla.t5, x = "treatment", y = "SLA", 
                 desc_stat = "mean_se", color = "black",
                 size=0.5,width = 2.5) +
  ylab('SLA') +
  xlab('Nutrient treatment') + 
  scale_y_continuous(breaks =seq(400,1600,200)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "SLA",
  width = 10,
  height = 10,
  noRStudioGD = TRUEf
)
SLA


#######################################################
#######################################################
#######################################################
### --- Goal 3.) Effects on Defensive Chemistry --- ###
#######################################################
#######################################################
#######################################################
# QUESTION: Does soil N and P availability affect a P. biflora defensive trait (leaf cyanide concentration)?
# HYPOTHESIS: cyanide concentration will increase linearly with N
# PLOTS:errorplots - average values plotted with standard error bars on both sides

### --- subset just final data collection 
t5 <- data %>% filter(timestamp == 5)

### --- This Rdata file has data objects for making all richness and diversity figures --- ###
dats=load(file='data/metabolomics/NPclass_groups.Rdata')
dats

### //////////////////////////////////////////////////////////////////////////////// ###
### /// FIGURE 3. Effects on P.biflora Chemical Traits - Defensive & Metabolomic /// ###
### //////////////////////////////////////////////////////////////////////////////// ###

# ------------------------------- #
# --- FIG. 3A: leaf N content --- #
# ------------------------------- #
N.t5<-t5[!is.na(t5$leafN),]
N.t5=N.t5[-86,] # remove the outlier to conform to normality (ESI processing error)
N.t5$treatment <- factor(N.t5$treatment,
                         levels = c("LNLP", "LNMP", "LNHP",
                                    "MNLP", "MNMP", "MNHP",
                                    "HNLP", "HNMP", "HNHP"))

N<-ggerrorplot(N.t5, x = "treatment", y = "leafN", 
               desc_stat = "mean_se", color = "black",
               size=0.5,width = 2.5) +
  ylab('Leaf N (%)') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "leaf N",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
N


# ------------------------------------------------- #
# --- FIG. 3B: N-containing metabolite richness --- #
# ------------------------------------------------- #
n_rich<-n_rich[!is.na(n_rich$treatment),]

n_rich$treatment <- factor(n_rich$treatment,
                            levels = c("LNLP", "LNMP", "LNHP",
                                       "MNLP", "MNMP", "MNHP",
                                       "HNLP", "HNMP", "HNHP"))

Ncontain<-ggerrorplot(n_rich, x = "treatment", y = "n", 
                      desc_stat = "mean_se", color = "black",
                      size=0.5,width = 2.5) +
  ylab('Nitrogen-containing metabolite richness') +
  xlab('Nutrient treatment') + 
  # ylim(150,200) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "N-containing richness",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
Ncontain


# --------------------------------- #
# --- FIG. 3C: CG concentration --- #
# --------------------------------- #

CG.t5<-t5[!is.na(t5$ug_mg),]
CG.t5=CG.t5[-93,] # remove outlier to conform to norm assumption (it doesn't change the LME or posthoc results)

CG.t5$treatment <- factor(CG.t5$treatment,
                         levels = c("LNLP", "LNMP", "LNHP",
                                    "MNLP", "MNMP", "MNHP",
                                    "HNLP", "HNMP", "HNHP"))

CGs<-ggerrorplot(CG.t5, x = "treatment", y = "ug_mg", 
               desc_stat = "mean_se", color = "black",
               size=0.5,width = 2.5) +
  ylab('Cyanogenic glycosides [ug/mg]') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "cyanide concentration",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
CGs


# ---------------------------------------------- #
# --- FIG. 3D: Cyanogenic glycoside richness --- #
# ---------------------------------------------- #
CG_rich<-CG_rich[!is.na(CG_rich$treatment),]

CG_rich$treatment <- factor(CG_rich$treatment,
                            levels = c("LNLP", "LNMP", "LNHP",
                                       "MNLP", "MNMP", "MNHP",
                                       "HNLP", "HNMP", "HNHP"))

CGspread<-ggerrorplot(CG_rich, x = "treatment", y = "n", 
                      desc_stat = "mean_se", color = "black",
                      size=0.5,width = 2.5) +
  ylab('Cyanogenic glycoside richness') +
  xlab('Nutrient treatment') + 
  ylim(0,4) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "cyanide richness",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
CGspread


# ------------------------------------------- #
# --- FIG. 3E: Total Metabolomic Richness --- #
# ------------------------------------------- #
met_rich<-met_rich[!is.na(met_rich$treatment),]

met_rich$treatment <- factor(met_rich$treatment,
                            levels = c("LNLP", "LNMP", "LNHP",
                                       "MNLP", "MNMP", "MNHP",
                                       "HNLP", "HNMP", "HNHP"))

metrico<-ggerrorplot(met_rich, x = "treatment", y = "n", 
                      desc_stat = "mean_se", color = "black",
                      size=0.5,width = 2.5) +
  ylab('Metabolomic richness') +
  xlab('Nutrient treatment') + 
  # ylim() +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "total metabolomic richness",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
metrico


# ---------------------------------------------------- #
# --- FIG. 3F: Total Biochemical Pathway Diversity --- #
# ---------------------------------------------------- #
path_div<-path_div[!is.na(path_div$treatment),]

path_div$treatment <- factor(path_div$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

bioDiv<-ggerrorplot(path_div, x = "treatment", y = "val", 
                     desc_stat = "mean_se", color = "black",
                     size=0.5,width = 2.5) +
  ylab('Biochemical pathway diversity') +
  xlab('Nutrient treatment') + 
  ylim(1.7,1.9) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "biochemical pathway diversity",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
bioDiv


### //////////////////////////////////////////////////////////////////////// ###
### /// FIGURE 4. NMDS Ordination Figure of Total Metabolomic Similarity /// ###
### //////////////////////////////////////////////////////////////////////// ###
### load the chemical similarity data 
dats=load(file="data/meabolomics/CSCS_biflora_nutrient-secondary_20221026.RData")
dats
# [1] "braydist"     "braydist.bin" "pd.class"     "class"        "date"         "meta" 

### Remove species pools from matrix
dim(braydist.bin) # [1] 147 147
braydist.bin.samps = braydist.bin[,-grep("_pool",names(braydist.bin))] # remove columns with pools
row.names(braydist.bin.samps)
braydist.bin.samps = braydist.bin.samps[-c(48,126,127,133,136,138,140,142,143),] # remove rows with pools 
dim(braydist.bin.samps) # [1] 138 138

### calculate NMDS using MASS
#nmds.samps = isoMDS(as.dist(braydist.bin.samps), k=5)
#nmds.samps.points = as.data.frame(nmds.samps$points)
#nmds.ready<-merge(nmds.samps.points,meta[,c(3:6)],by="ATTRIBUTE_SampleNumber")
#grp.h <- nmds.ready[nmds.ready$ATTRIBUTE_Nitrogen == "high", ][chull(nmds.ready[nmds.ready$ATTRIBUTE_Nitrogen == 
#                                                                                  "high", c("V1", "V2")]), ]  # hull values for grp A
#grp.m <- nmds.ready[nmds.ready$ATTRIBUTE_Nitrogen == "medium", ][chull(nmds.ready[nmds.ready$ATTRIBUTE_Nitrogen == 
#                                                                                    "medium", c("VS", "V2")]), ]  # hull values for grp B
#grp.l <- nmds.ready[nmds.ready$ATTRIBUTE_Nitrogen == "low", ][chull(nmds.ready[nmds.ready$ATTRIBUTE_Nitrogen ==                                                                              "low", c("V1", "V2")]), ]  # hull values for grp B
#hull.data <- rbind(grp.h, grp.m, grp.l)  #combine grp.a and grp.b
#nmds.samps.points$ATTRIBUTE_SampleNumber=row.names(nmds.samps.points)

### calculate NMDS axes using vegan
library(vegan)
nmds.samps2 = metaMDS(as.dist(braydist.bin.samps), k=5)
# TOTAL Metabolome Stress:     0.08423414 
# SECONDARY METABOLITE Stress:   0.08331218 
# PRIMARY METABOLITE Stress:  0.09407974 
nmds.samps.points2 = as.data.frame(nmds.samps2$points)

### add column with sample names to NMDS points df for merging with meta data file
nmds.samps.points2$ATTRIBUTE_SampleNumber=row.names(nmds.samps.points2)

### load metadata so that we can link sample numbers with their treatments 
#meta = read.delim("data/metabolomics/MZmine_output/G039_sampspools_metadata.txt", header = T, sep = "	", comment.char = "")

### merge with treatment categories from meta data file
nmds.ready2<-merge(nmds.samps.points2,meta[,c(3:6)],by="ATTRIBUTE_SampleNumber")
head(nmds.ready2)

### make another object grouped by treatments for the polygon hulls 
grp.h2 <- nmds.ready2[nmds.ready2$ATTRIBUTE_Nitrogen == "high", ][chull(nmds.ready2[nmds.ready2$ATTRIBUTE_Nitrogen == 
                                                                                      "high", c("MDS1", "MDS2")]), ]  # hull values for grp A
grp.m2 <- nmds.ready2[nmds.ready2$ATTRIBUTE_Nitrogen == "medium", ][chull(nmds.ready2[nmds.ready2$ATTRIBUTE_Nitrogen == 
                                                                                        "medium", c("MDS1", "MDS2")]), ]  # hull values for grp B
grp.l2 <- nmds.ready2[nmds.ready2$ATTRIBUTE_Nitrogen == "low", ][chull(nmds.ready2[nmds.ready2$ATTRIBUTE_Nitrogen == 
                                                                                     "low", c("MDS1", "MDS2")]), ]  # hull values for grp B
hull.data2 <- rbind(grp.h2, grp.m2, grp.l2)  #combine grp.a and grp.b
hull.data2

### make a color scheme
library(viridis)
colors = viridis(10) # [1] "#440154FF" "#22A884FF" "#B4DE2CFF"

### plot it
project<-ggplot() + 
  geom_point(data=nmds.ready2,aes(x=MDS1,y=MDS2,colour=ATTRIBUTE_Nitrogen),size=3) + # add the point markers
  geom_polygon(data=hull.data2,aes(x=MDS1,y=MDS2,fill=ATTRIBUTE_Nitrogen,group=ATTRIBUTE_Nitrogen),alpha=0.2) + # add the convex hulls
  scale_colour_manual(values=c("high" = "#440154FF", "medium" = "#22A884FF", "low" = "#B4DE2CFF")) +
  scale_fill_manual(values=c("high" = "#440154FF", "medium" = "#22A884FF", "low" = "#B4DE2CFF")) +
  labs(colour="Nitrogen treatment",
       fill="Nitrogen treatment",
       x="NMDS axis 1",
       y="NMDS axis 2") +
  coord_equal() +
  theme_bw() +
  theme(legend.title = element_text(size=16,face="bold"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        legend.text = element_text(size=16,color = "black"),
        axis.text.x = element_text(size=14,color = "black"),
        axis.text.y = element_text(size=14,color = "black"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = -15, unit = "pt")
  )
dev.new(
  title = "metabolome NMDS",
  width = 10,
  height = 7,
  noRStudioGD = TRUE
)
project



####################################
####################################
### ---------------------------- ###
### --- Supplemental Figures --- ###
### ---------------------------- ###
####################################
####################################

# ------------------------------------------------------- #
# --- FIG. S6: Plot counts of samples with unique CGs --- #
# ------------------------------------------------------- #
dats=load(file='data/metabolomics/NPclass_groups.Rdata')
dats

CG_count
# subset each unique CG and order the factor levels appropriately 
ali412<-subset(CG_count, row.ID == 412)
ali412$Nlevel <- factor(ali412$Nlevel,
                        levels = c("low", "medium", "high"))

ali1322<-subset(CG_count, row.ID == 1322)
ali1322$Nlevel <- factor(ali1322$Nlevel,
                         levels = c("low", "medium", "high"))

ali365<-subset(CG_count, row.ID == 365)
ali365$Nlevel <- factor(ali365$Nlevel,
                        levels = c("low", "medium", "high"))

aro333<-subset(CG_count, row.ID == 333)
aro333$Nlevel <- factor(aro333$Nlevel,
                        levels = c("low", "medium", "high"))

####  Four bargraphs
# compound 333
a412<-ggplot(ali412, aes(x=Nlevel, y=n)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           fill = "grey70",
           position=position_dodge(),
           lwd=0.75) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  ylim(0,50) + 
  #geom_text(size=5,
  #          x = 1, y = 51,label= "aliphatic CG") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=2),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face='bold'),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
a412

# compound 333
b1322<-ggplot(ali1322, aes(x=Nlevel, y=n)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           fill = "grey70",
           position=position_dodge(),
           lwd=0.75) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  ylim(0,50) + 
  #geom_text(size=5,
  #          x = 1, y = 51,label= "aliphatic CG") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=2),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face='bold'),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
b1322

# compound 333
c365<-ggplot(ali365, aes(x=Nlevel, y=n)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           fill = "grey70",
           position=position_dodge(),
           lwd=0.75) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  ylim(0,50) + 
  #geom_text(size=5,
  #          x = 1, y = 51,label= "aliphatic CG") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=2),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face='bold'),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
c365

# compound 333
d333<-ggplot(aro333, aes(x=Nlevel, y=n)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           fill = "grey70",
           position=position_dodge(),
           lwd=0.75) +
  ggtitle("") +
  xlab("") + 
  ylab("") +
  ylim(0,50) + 
  #geom_text(size=5,
  #          x = 1, y = 51,label= "aromatic CG") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 20, face = "bold", colour = "black"),
        legend.title = element_text(size = 18, face = "bold", colour = "black"),
        legend.text = element_text(size = 18, colour = "black"),
        legend.background = element_rect(size=0.25, linetype="solid", 
                                         colour ="black"),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size=18),
        axis.text.x = element_text(size=18),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=18, face='bold',vjust=2),
        axis.title.x = element_text(margin = margin(t = 20),size=18,face='bold'),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
d333

### --- make survival over time  plots into multipanel figure
comps=ggarrange(a412,b1322,
               c365,d333,
               labels = c("(a)", "(b)",
                          "(c)","(d)"),
               font.label = list(size = 18),
               nrow=2,ncol=2)
# dev plot
dev.new(
  title = "comp.counts",
  width = 8.5,
  height = 11,
  noRStudioGD = TRUE
  )
comps

annotate_figure(comps,
                bottom = text_grob("Nitrogen level", vjust=-1.0,hjust=0.25, size = 18,face='bold'),
                left = text_grob("Samples with this metabolite", vjust=0.25,hjust=0.3, size = 18,rot=90,face='bold'))



# ---------------------------------------------------------------------- #
# --- FIG. S7: Richness & pathway diversity by primary and secondary --- #
# ---------------------------------------------------------------------- #
### A: primary richness
prim_rich<-prim_rich[!is.na(prim_rich$treatment),]

prim_rich$treatment <- factor(prim_rich$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

primrico<-ggerrorplot(prim_rich, x = "treatment", y = "n", 
                     desc_stat = "mean_se", color = "black",
                     size=0.5,width = 2.5) +
  ylab('Primary richness') +
  xlab('') + 
  # ylim(100,280) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,face='bold'),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )
primrico

### B: secondary richness
second_rich<-second_rich[!is.na(second_rich$treatment),]
second_rich$treatment <- factor(second_rich$treatment,
                              levels = c("LNLP", "LNMP", "LNHP",
                                         "MNLP", "MNMP", "MNHP",
                                         "HNLP", "HNMP", "HNHP"))

segunrico<-ggerrorplot(second_rich, x = "treatment", y = "n", 
                      desc_stat = "mean_se", color = "black",
                      size=0.5,width = 2.5) +
  ylab('Secondary richness') +
  xlab('') + 
  # ylim(100,280) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,face='bold'),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )
segunrico

### C: primary biochemical pathway diversity 
prim_div<-prim_div[!is.na(prim_div$treatment),]
prim_div$treatment <- factor(prim_div$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

primDiv<-ggerrorplot(prim_div, x = "treatment", y = "val", 
                    desc_stat = "mean_se", color = "black",
                    size=0.5,width = 2.5) +
  ylab('Primary diversity') +
  xlab('') + 
  #ylim(0.25,1.7) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,face='bold'),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )
primDiv

### D: secondary biochemical pathway diversity 
second_div<-second_div[!is.na(second_div$treatment),]
second_div$treatment <- factor(second_div$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

segunDiv<-ggerrorplot(second_div, x = "treatment", y = "val", 
                     desc_stat = "mean_se", color = "black",
                     size=0.5,width = 2.5) +
  ylab('Secondary diversity') +
  xlab('') + 
  # ylim(0.25,1.7) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,face="bold"),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )
segunDiv

### --- assemble the plots into a multipanel figure
breaks=ggarrange(primrico,segunrico,
                 primDiv,segunDiv,
                 labels = c("(a)", "(b)",
                            "(c)","(d)"),
                 font.label = list(size = 18),
                 nrow=2,ncol=2)
dev.new(
  title = "richness and div broken down",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
breaks

annotate_figure(breaks,
                bottom = text_grob("Nutrient treatment", vjust=-1.0,hjust=0.25, size = 18,face='bold'))








###########################################
###########################################
### ----------------------------------- ###
### ---  Figures NOT Included in MS --- ###
### ----------------------------------- ###
###########################################
###########################################

# ---------------------- #
# --- leaf C content --- #
# ---------------------- #
C.t5<-t5[!is.na(t5$leafC),]
C.t5=C.t5[-86,] # N.t5=N.t5[-86,] # remove the outlier to conform to normality (ESI processing error)

C.t5$treatment <- factor(C.t5$treatment,
                         levels = c("LNLP", "LNMP", "LNHP",
                                    "MNLP", "MNMP", "MNHP",
                                    "HNLP", "HNMP", "HNHP"))
C<-ggerrorplot(CN.t5, x = "treatment", y = "leafC", 
               desc_stat = "mean_se", color = "black",
               size=0.5,width = 2.5) +
  ylab('Leaf C (%)') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "carbon",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
C


# --------------- #
# --- soil pH --- #
# --------------- #
pH.t5<-t5[!is.na(t5$pH),]
pH.t5$treatment <- factor(pH.t5$treatment,
                          levels = c("LNLP", "LNMP", "LNHP",
                                     "MNLP", "MNMP", "MNHP",
                                     "HNLP", "HNMP", "HNHP"))

pH<-ggerrorplot(t5, x = "treatment", y = "pH", 
                desc_stat = "mean_se", color = "black",
                size=0.5,width = 2.5) +
  ylab('Soil pH') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(400,1600,200)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "pH",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
pH


# -------------------- #
# --- leafCN_ratio --- #
# -------------------- #
CN.t5<-t5[!is.na(t5$leafCN_ratio),]
CN.t5=CN.t5[-86,] # N.t5=N.t5[-86,] # remove the outlier to conform to normality (ESI processing error)

CN.t5$treatment <- factor(CN.t5$treatment,
                          levels = c("LNLP", "LNMP", "LNHP",
                                     "MNLP", "MNMP", "MNHP",
                                     "HNLP", "HNMP", "HNHP"))

CN<-ggerrorplot(CN.t5, x = "treatment", y = "leafCN_ratio", 
                desc_stat = "mean_se", color = "black",
                size=0.5,width = 2.5) +
  ylab('Carbon/Nitrogen ratio') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "root:shoot",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
CN


# ------------------------ #
# --- Internode length --- #
# ------------------------ #
node.t5<-t5[!is.na(t5$internode_length),]
node.t5$treatment <- factor(node.t5$treatment,
                            levels = c("LNLP", "LNMP", "LNHP",
                                       "MNLP", "MNMP", "MNHP",
                                       "HNLP", "HNMP", "HNHP"))

node<-ggerrorplot(node.t5, x = "treatment", y = "internode_length", 
                  desc_stat = "mean_se", color = "black",
                  size=0.5,width = 2.5) +
  ylab('Internode length (cm)') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "internode",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
node


# ----------------- #
# --- Leaf Area --- #
# ----------------- #
Adata <- read_xlsx("SLA_21area.xlsx")
At5 <- Adata %>% filter(timestamp == 5)
At5$leaf_area<-as.numeric(as.character(At5$leaf_area))
At5<-subset(At5, !is.na(leaf_area))

At5$treatment <- factor(At5$treatment,
                        levels = c("LNLP", "LNMP", "LNHP",
                                   "MNLP", "MNMP", "MNHP",
                                   "HNLP", "HNMP", "HNHP"))

area<-ggerrorplot(At5, x = "treatment", y = "leaf_area", 
                  desc_stat = "mean_se", color = "black",
                  size=0.5,width = 2.5) +
  ylab(expression(paste("Leaf area (cm"^"2",")"))) +
  xlab('Nutrient treatment') + 
  scale_y_continuous(breaks =seq(400,1600,200)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "leaf area",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
area


# ------------------ #
# --- mass_below --- #
# ------------------ #
below.t5<-t5[!is.na(t5$mass_below),]
below.t5$treatment <- factor(below.t5$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

below<-ggerrorplot(below.t5, x = "treatment", y = "mass_below", 
                   desc_stat = "mean_se", color = "black",
                   size=0.5,width = 2.5) +
  ylab('Belowground biomass (g)') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "below biomass",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
below


# ------------------ #
# --- mass_above --- #
# ------------------ #
above.t5<-t5[!is.na(t5$mass_above),]
above.t5$treatment <- factor(above.t5$treatment,
                             levels = c("LNLP", "LNMP", "LNHP",
                                        "MNLP", "MNMP", "MNHP",
                                        "HNLP", "HNMP", "HNHP"))

above<-ggerrorplot(below.t5, x = "treatment", y = "mass_below", 
                   desc_stat = "mean_se", color = "black",
                   size=0.5,width = 2.5) +
  ylab('Aboveground biomass (g)') +
  xlab('Nutrient treatment') + 
  # scale_y_continuous(breaks =seq(0,2.0,0.25)) +
  geom_vline(xintercept = 3.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  geom_vline(xintercept = 6.5, linetype="longdash", 
             color = "grey50", size=0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 30,hjust=1,size=12,color = "black"),
        axis.text.y = element_text(size=12,color = "black"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16),
        plot.margin = margin(t = 10, r = 10, b = -10, l =  -10, unit = "pt")
  )

dev.new(
  title = "above biomass",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
above



###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################