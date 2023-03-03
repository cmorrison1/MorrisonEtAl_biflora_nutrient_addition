
library(phytools)
library(picante)
library(MASS)
library(viridis)
library(RColorBrewer)
library(vegan)
library(caper)

getwd()
setwd('~/Desktop/biflora.nutrient.addition')

gnps.code.feat = "2014c4a5" #demo_10k Tyson Reearch Center hickories (Carya)
gnps.code.npclass = "2ff7c4d9"

feat = read.csv("data/metabolomics/FBMN_output/ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-2014c4a5-download_qza_table_data/quantification_table/quantification_table-00000.csv")

dim(feat)
# [1] 2083  154
head(feat)

meta = read.delim("data/metabolomics/MZmine_output/G039_sampspools_metadata.txt", header = T, sep = "	", comment.char = "")
dim(meta)
# [1] 150 9

tree = read.tree("data/metabolomics/SiriusQemistreeNPClass_output/G039_samps_20k_qemistree.nwk")
# Phylogenetic tree with 973 tips and 972 internal nodes.

siri = read.table("data/metabolomics/SiriusQemistreeNPClass_output/G039_samps_20k_classified_feature_data.tsv", header = T, sep = "	", comment.char = "", quote = "")
names(siri)[1] = "id"
dim(siri)
# [1] 973  15

npclass = read.table(paste("data/metabolomics/SiriusQemistreeNPClass_output/ProteoSAFe-NPCLASSIFIER-", gnps.code.npclass, "-view_results/NPCLASSIFIER-", gnps.code.npclass, "-view_results-main.tsv", sep = ""), header = T, sep = "\t", comment.char = "")
#npclass = read.table("data/metabolomics/SiriusQemistreeNPClass_output/ProteoSAFe-NPCLASSIFIER-2ff7c4d9-view_results/NPCLASSIFIER-2ff7c4d9-view_results-main.tsv", header = T, sep = "\t", comment.char = "")
dim(npclass)
#[1] 785   5
head(npclass)

# smiles             class_results superclass_results
# 1                                      COC1=CC(=CC2=C1C(=O)C=C(O2)C3=CC=CC=C3O)O                  Flavones         Flavonoids
# 2                                   CC(C)C(C(=O)OCC(CO)OCN1C=NC2=C1N=C(NC2=O)N)N       pteridine alkaloids    Pseudoalkaloids
# 3                         C1=C(C=C2C(=C1O)C(=O)C(=C(O2)C3=CC(=C(C(=C3O)O)O)O)O)O                 Flavonols         Flavonoids
# 4                                  C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O                 Flavonols         Flavonoids
# 5                                             CCCCCC=CCC=CCC=CCC1C(O1)CCCC(=O)OC Epoxyeicosatrienoic acids        Eicosanoids
# 6 CC(=O)OC1=C(C=C(C=C1)C2=C(C(=O)C3=C(O2)C=C(C=C3OC(=O)C)OC(=O)C)OC(=O)C)OC(=O)C                 Flavonols         Flavonoids
# pathway_results isglycoside
# 1 Shikimates and Phenylpropanoids       False
# 2                       Alkaloids       False
# 3 Shikimates and Phenylpropanoids       False
# 4 Shikimates and Phenylpropanoids       False
# 5                     Fatty acids       False
# 6 Shikimates and Phenylpropanoids       False


siri$class_results = NA
siri$superclass_results = NA
siri$pathway_results = NA
siri$isglycoside = NA

for(i in 1:nrow(siri)){
  smiles_i = siri$smiles[i]
  if(smiles_i %in% npclass$smiles){
    siri$class_results[i] = npclass$class_results[which(npclass$smiles == smiles_i)][1]
    siri$superclass_results[i] = npclass$superclass_results[which(npclass$smiles == smiles_i)][1]
    siri$pathway_results[i] = npclass$pathway_results[which(npclass$smiles == smiles_i)][1]
    siri$isglycoside[i] = npclass$isglycoside[which(npclass$smiles == smiles_i)][1]
  }
}

dim(siri)
# [1] 973  19

head(siri)

# #juglone
# grep("O=C\2c1c(c(O)ccc1)C(=O)/C=C/2", siri$smiles)
# #Polyketides: Naphthalenes: Naphthoquinones
# siri[which(siri$class_results == "Naphthoquinones"),]
# COC1=C2C=CC=C(C2=C(C=C1)OC)OC ## this is showing up in TRC Carya and Juglans (methylated juglone)

#siri$id[which(siri$smiles == "COC1=C2C=CC=C(C2=C(C=C1)OC)OC")]
# [1] "25efb88e3fe244a70a23437e88bab80b"

#siri[which(siri$smiles == "COC1=C2C=CC=C(C2=C(C=C1)OC)OC"),]



siri$X.featureID[i]
# [1] "700,596"

#for(i in 1:nrow(siri)){
#  if(length(strsplit(siri$X.featureID[i], ",")[[1]]) > 1){
#    featureID = strsplit(siri$X.featureID[i], ",")[[1]]
#    for(j in 1:length(featureID)){
#      siri$parent_mass[i] = feat$row.m.z[which(feat$row.ID == featureID[j])]
#      siri$retention_time[i] = feat$row.retention.time[which(feat$row.ID == featureID[j])]
#    }
#  }
#  if(length(strsplit(siri$X.featureID[i], ",")[[1]]) == 1){
#    featureID = siri$X.featureID[i]
#    siri$parent_mass[i] = feat$row.m.z[which(feat$row.ID == featureID)]
#    siri$retention_time[i] = feat$row.retention.time[which(feat$row.ID == featureID)]
#  }
#  
#}
for(i in 1:nrow(siri)){
  featureID = siri$X.featureID[i]
  siri$parent_mass[i] = feat$row.m.z[which(feat$row.ID == featureID)]
  siri$retention_time[i] = feat$row.retention.time[which(feat$row.ID == featureID)]
}



nsamps = ncol(feat)-4
nsamps
# [1] 150


heat = as.data.frame(matrix(0,nrow = nrow(siri), ncol = (nsamps)))
row.names(heat) = siri$id
names(heat) = names(feat)[4:(ncol(feat)-1)]

dim(heat)
# [1] 973 150

### if merging multiple
multips = grep(",", siri$table_number)
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

dim(heat)
# [1] 973 150

# length(which(meta$ATTRIBUTE_SampleType == "BLANK"))

# grep("LANK", names(heat))
# which(meta$ATTRIBUTE_SampleType == "BLANK")
# meta$filename[which(meta$ATTRIBUTE_SampleType == "BLANK")]

grep("lank", names(heat))
# [1] 1   2   114

# blanklist = c(meta$filename[which(meta$ATTRIBUTE_SampleType == "BLANK")])
# blankspots = rep(0,length(which(meta$ATTRIBUTE_SampleType == "BLANK")))
# for(i in 1:length(blanklist)){
# blankspots[i] = grep(blanklist[i],names(heat))
# }
# blankspots
# heat.real = heat[which(rowSums(heat[,blankspots]) == 0),]


heat.real = heat[which(rowSums(heat[,grep("lank", names(heat))]) == 0),]
#View(heat.real)
dim(heat.real)
# [1] 735 150

colSums(heat.real)


# colSums(heat.real)[blankspots]
# heat.real = heat.real[,-blankspots]

heat.real = heat.real[,-grep("lank", names(heat))]

dim(heat.real)
# [1] 735 147

siri.real = siri[which(siri$id %in% row.names(heat.real)),]

dim(siri.real)
# [1] 735  19


levels(as.factor(siri.real$pathway_results))
# [1] ""                                                        
# [2] "Alkaloids"                                               
# [3] "Alkaloids,Amino acids and Peptides"                      
# [4] "Alkaloids,Shikimates and Phenylpropanoids"               
# [5] "Amino acids and Peptides"                                
# [6] "Amino acids and Peptides,Fatty acids"                    
# [7] "Amino acids and Peptides,Shikimates and Phenylpropanoids"
# [8] "Carbohydrates"                                           
# [9] "Fatty acids"                                             
# [10] "Polyketides"                                             
# [11] "Shikimates and Phenylpropanoids"                         
# [12] "Terpenoids" 


siri.real$custom = NA
siri.real$custom[which(siri.real$pathway_results == "Alkaloids")] = "Alkaloids"
siri.real$custom[which(siri.real$pathway_results == "Alkaloids,Amino acids and Peptides")] = "Alkaloids"
siri.real$custom[which(siri.real$pathway_results == "Alkaloids,Carbohydrates")] = "Alkaloids"
siri.real$custom[which(siri.real$pathway_results == "Alkaloids,Fatty acids")] = "Alkaloids"

siri.real$custom[which(siri.real$pathway_results == "Polyketides")] = "Polyketides"
siri.real$custom[which(siri.real$pathway_results == "Amino acids and Peptides,Polyketides")] = "Polyketides"

siri.real$custom[which(siri.real$pathway_results == "Shikimates and Phenylpropanoids")] = "Shikimates and Phenylpropanoids"
siri.real$custom[which(siri.real$pathway_results == "Amino acids and Peptides,Shikimates and Phenylpropanoids")] = "Shikimates and Phenylpropanoids"

siri.real$custom[which(siri.real$pathway_results == "Terpenoids")] = "Terpenoids"

siri.real$custom[which(siri.real$pathway_results == "Alkaloids,Polyketides")] = "Multiple"
siri.real$custom[which(siri.real$pathway_results == "Alkaloids,Shikimates and Phenylpropanoids")] = "Multiple"
siri.real$custom[which(siri.real$pathway_results == "Alkaloids,Terpenoids")] = "Multiple"
siri.real$custom[which(siri.real$pathway_results == "Polyketides,Shikimates and Phenylpropanoids")] = "Multiple"
siri.real$custom[which(siri.real$pathway_results == "Polyketides,Terpenoids")] = "Multiple"
siri.real$custom[which(siri.real$pathway_results == "Shikimates and Phenylpropanoids,Terpenoids")] = "Multiple"

siri.real$custom[which(siri.real$pathway_results == "Amino acids and Peptides")] = "Amino acids and Peptides"

siri.real$custom[which(siri.real$pathway_results == "Carbohydrates")] = "Carbohydrates"

siri.real$custom[which(siri.real$pathway_results == "Fatty acids")] = "Fatty acids"


head(siri.real)
head(heat.real)

demo.metab = cbind(siri.real, heat.real)
head(demo.metab)

# save(demo.metab, heat.real, siri.real, heat, meta, feat, tree, siri, npclass, file = "demo_carya_metab_20220822.RData")


demo.metab[which(demo.metab$smiles == "COC1=C2C=CC=C(C2=C(C=C1)OC)OC"),]



#### put prepwork for iTOL figure here, followed by calculation of CSCS and Bray-Curtis (and presence/absence similarity) for each chemical class
## note: use nodecolormat to make file to color qemistree phylogeny branches based on chemical class

## then calculate relative abundance of every tree species in every plot, use those numbers in the multibarplot, use tip.label to refer to tips for iTOL

# LEGEND_TITLE,Chemical Superclasses
# LEGEND_SHAPES,1,1,1,1,1
# LEGEND_LABELS,Alkaloids,Benzenoids,Heterocyclics,Lipids,Organic Acids
# LEGEND_COLORS,#00ff00,#0000ff,#ff0000,#ffff00,#00ffff

# 9afb472f40aaa5497b6ac48a8fe6db5c,branch,#00ff00,normal,1
# 380245445bf13e84a18fdae75c2804a8,branch,#00ffff,normal,1
# c0194ee4b6de5149f7c50365f8a9eab9,branch,#00ffff,normal,1
# 4ff24827c34deb8787732b3597007fd6,branch,#00ffff,normal,1
# 54dff4bfcdd92d66a4d037604afecf32,branch,#00ffff,normal,1




library(RColorBrewer)
cols.new = brewer.pal(8, "Set3")

cols.new
# [1] "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5"


siri.itol = siri.real
dim(siri.itol)
# [1] 735  20

levels(as.factor(siri.itol$custom))
# [1] "Alkaloids"                       "Amino acids and Peptides"        "Carbohydrates"                   "Fatty acids"                    
# [5] "Polyketides"                     "Shikimates and Phenylpropanoids" "Terpenoids"  


## In the BCI 50-ha forest plot dataset:
# [1] "Alkaloids"                       "Amino acids and Peptides"        "Carbohydrates"                   "Fatty acids"                    
# [5] "Multiple"                        "Polyketides"                     "Shikimates and Phenylpropanoids" "Terpenoids"  

length(which(siri.itol$custom == "NA"))
# [1] 0



comm.mat = as.data.frame(matrix(1, nrow=2, ncol = nrow(siri.itol)))
names(comm.mat) = siri.itol$id
qtree.itol =  prune.sample(samp = comm.mat, phylo = tree)

qtree.itol
# Phylogenetic tree with 735 tips and 734 internal nodes.

# Tip labels:
# 11e48d5cc6cecabf93148a6b1e828d71, fe6191242984742887c1dafc7258768b, 3482462162c3f7619ae60e8da89ae25d, 824e19749aa74c4bc65a19e9e9c537c6, 37941f0498adf6b33dd835fe10620191, 21476e26b1cd807ee83f6c19f2f551da, ...

# Rooted; includes branch lengths.

heat = heat.real
dim(heat)
# [1] 735 147

heat.itol = heat
# heat.itol = heat.itol[,-grep("TRP",names(heat.itol))]
heat.itol = heat.itol[qtree.itol$tip.label,]

rownames(siri.itol) = siri.itol$id
siri.itol = siri.itol[qtree.itol$tip.label,]

dim(siri.itol)
# [1] 735  20

itol.branchcolors = as.data.frame(matrix(NA,nrow = nrow(heat.itol), ncol = 5))
itol.branchcolors[,1] = siri.itol$id
itol.branchcolors[,2] = "branch"
itol.branchcolors[,4] = "normal"
itol.branchcolors[,5] = 1

head(itol.branchcolors)
#                                V1     V2 V3     V4 V5
#1 09502dfff5a9821f783515845fb4203b branch NA normal  1
#2 5a4a48a0ed76f2b2c65e8de1690416f1 branch NA normal  1
#3 c5dad8482e299229ad104bc35ef2dc06 branch NA normal  1
#4 5fccb31efda4962846b6693a07a76ec5 branch NA normal  1
#5 bffc55d3cb6b06dfd85a033b6ebe8bbb branch NA normal  1
#6 2bcf91a01e154e5293406195295e5e7c branch NA normal  1

row.names(itol.branchcolors) = itol.branchcolors$V1

itol.branchcolors[which(siri.itol$custom == "Alkaloids"),3] = cols.new[1]
itol.branchcolors[which(siri.itol$custom == "Amino acids and Peptides"),3] = cols.new[2]
itol.branchcolors[which(siri.itol$custom == "Carbohydrates"),3] = cols.new[3]
itol.branchcolors[which(siri.itol$custom == "Fatty acids"),3] = cols.new[4]
itol.branchcolors[which(siri.itol$custom == "Polyketides"),3] = cols.new[5]
itol.branchcolors[which(siri.itol$custom == "Shikimates and Phenylpropanoids"),3] = cols.new[6]
itol.branchcolors[which(siri.itol$custom == "Terpenoids"),3] = cols.new[7]
itol.branchcolors[which(siri.itol$custom == "Multiple"),3] = cols.new[8]

cols.new
# [1] "#8DD3C7" "#FFFFB3" "#BEBADA" "#FB8072" "#80B1D3" "#FDB462" "#B3DE69" "#FCCDE5"


head(itol.branchcolors)
# V1     V2      V3     V4 V5
# 09502dfff5a9821f783515845fb4203b 09502dfff5a9821f783515845fb4203b branch    <NA> normal  1
# 5a4a48a0ed76f2b2c65e8de1690416f1 5a4a48a0ed76f2b2c65e8de1690416f1 branch #FFFFB3 normal  1
# c5dad8482e299229ad104bc35ef2dc06 c5dad8482e299229ad104bc35ef2dc06 branch #FFFFB3 normal  1
# 5fccb31efda4962846b6693a07a76ec5 5fccb31efda4962846b6693a07a76ec5 branch #BEBADA normal  1
# bffc55d3cb6b06dfd85a033b6ebe8bbb bffc55d3cb6b06dfd85a033b6ebe8bbb branch #BEBADA normal  1
# 2bcf91a01e154e5293406195295e5e7c 2bcf91a01e154e5293406195295e5e7c branch #BEBADA normal  1

dim(itol.branchcolors)
# [1] 735   5

itol.branchcolors = itol.branchcolors[qtree.itol$tip.label,]


qtree.itol
# Phylogenetic tree with 735 tips and 734 internal nodes.



head(itol.branchcolors)
# V1     V2      V3     V4 V5
# 09502dfff5a9821f783515845fb4203b 09502dfff5a9821f783515845fb4203b branch    <NA> normal  1
# 5a4a48a0ed76f2b2c65e8de1690416f1 5a4a48a0ed76f2b2c65e8de1690416f1 branch #FFFFB3 normal  1
# c5dad8482e299229ad104bc35ef2dc06 c5dad8482e299229ad104bc35ef2dc06 branch #FFFFB3 normal  1
# 5fccb31efda4962846b6693a07a76ec5 5fccb31efda4962846b6693a07a76ec5 branch #BEBADA normal  1
# bffc55d3cb6b06dfd85a033b6ebe8bbb bffc55d3cb6b06dfd85a033b6ebe8bbb branch #BEBADA normal  1
# 2bcf91a01e154e5293406195295e5e7c 2bcf91a01e154e5293406195295e5e7c branch #BEBADA normal  1

# itol.branchcolors = itol.branchcolors[-which(is.na(itol.branchcolors[,3])),]

dim(itol.branchcolors)
# [1] 735   5


# write.csv(itol.branchcolors, file = "demo_carya_10k_NPClassifier_iTOLbranchcolors_table_20220822.csv", quote = FALSE, row.names = FALSE)

# write.tree(qtree.itol, file = "demo_carya_10k_qtree_pruned_20220822.tre")

# heat.itol = heat
# heat.itol = heat.itol[qtree.itol$tip.label,]
itol.barchart.sp = as.data.frame(matrix(NA,nrow = nrow(heat.itol), ncol = 2))
itol.barchart.sp[,1] = row.names(heat.itol)

for(i in 1:nrow(heat.itol))
{
  itol.barchart.sp[i,2] = length(which(heat.itol[i,] > 0))
}

itol.barchart.sp[1:10,]
# V1 V2
# 1  09502dfff5a9821f783515845fb4203b  30
# 2  5a4a48a0ed76f2b2c65e8de1690416f1 104
# 3  c5dad8482e299229ad104bc35ef2dc06 134
# 4  5fccb31efda4962846b6693a07a76ec5 111
# 5  bffc55d3cb6b06dfd85a033b6ebe8bbb  15
# 6  2bcf91a01e154e5293406195295e5e7c 131
# 7  8a69f1a666c67db755bad3631d1a7e73 139
# 8  bbd1afc1c43c5567da3e09847c6f4da5  37
# 9  fae3ecda58ef46661e5c905573ad15cf 146
# 10 4a30f587f58f10357438abdedadf6d22 137

# write.csv(itol.barchart.sp, file = "demo_carya_10k_NPClassifier_iTOLbarchartSpp_table_20220822.csv", quote = FALSE, row.names = FALSE)

### These determine metabolite class (or total for all metabolites) that will be summarized by CSCS, BC, etc.  
class = "total"
class = "secondary"
class = "primary"

class = "alkaloids"
class = "peptides"
class = "polyketides"
class = "shikimates"
class = "terpenoids"
class = "carbohydrates"
class = "fatty_acids"



date = 20230120
outfile = paste("CSCS_biflora_binCSCStest-", class, "_", date, ".RData",sep="")





# for(j in 1:ncol(heat.real)){
# if(length(which(heat.real[,j] > 0)) <2){
# cat("Sample", names(heat.real)[j], "has only 1 compound", "\n", sep = " ")
# }
# }

# "Sample TR00027_CARYTE.mzXML.Peak.area has only 1 compound" 

siri.itol$id[which(siri.itol$smiles == "COC1=C2C=CC=C(C2=C(C=C1)OC)OC")]
# [1] "25efb88e3fe244a70a23437e88bab80b"

siri.itol[which(siri.itol$smiles == "COC1=C2C=CC=C(C2=C(C=C1)OC)OC"),]



### select the chemical group you want to analyze
if(class == "total"){
  heat.class = heat.itol
}
if(class == "secondary"){
  heat.class = heat.itol[which(siri.itol$custom %in% c("Alkaloids","Amino acids and Peptides","Polyketides","Shikimates and Phenylpropanoids", "Terpenoids", "Multiple")),]
}
if(class == "primary"){
  heat.class = heat.itol[-which(siri.itol$custom %in% c("Alkaloids","Amino acids and Peptides","Polyketides","Shikimates and Phenylpropanoids", "Terpenoids", "Multiple")),]
}
if(class == "alkaloids"){
  heat.class = heat.itol[grep("Alkaloids", siri.itol$pathway_results),]
}
if(class == "peptides"){
  heat.class = heat.itol[grep("Peptides", siri.itol$pathway_results),]
}
if(class == "polyketides"){
  heat.class = heat.itol[grep("Polyketides", siri.itol$pathway_results),]
}
if(class == "shikimates"){
  heat.class = heat.itol[grep("Shikimates", siri.itol$pathway_results),]
}
if(class == "terpenoids"){
  heat.class = heat.itol[grep("Terpenoids", siri.itol$pathway_results),]
}
if(class == "fatty_acids"){
  heat.class = heat.itol[which(siri.itol$pathway_results == "Fatty acids"),]
}
if(class == "carbohydrates"){
  heat.class = heat.itol[which(siri.itol$pathway_results == "Carbohydrates"),]
}

# heat.class = heat.class[,-which(names(heat.class) == "TR00027_CARYTE.mzXML.Peak.area")] # ??????

sampsByCompounds = as.data.frame(t(as.data.frame(heat.class)))
names(sampsByCompounds) = row.names(heat.class)

comm.mat = as.data.frame(matrix(1, nrow=2, ncol = ncol(sampsByCompounds)))
names(comm.mat) = names(sampsByCompounds)
qtree.class = prune.sample(samp = comm.mat, phylo = qtree.itol)
qtree.class = qtree.itol
# maxdist = max(cophenetic(qtree.class))


#### CSCS loop starts here 
maxdist = pd.calc(cm = tree, tip.subset = c(tree$tip.label[1], tree$tip.label[length(tree$tip.label)]))[1]
pairwise.comps = cophenetic(qtree.class)
pairwise.comps = 1-(cophenetic(qtree.class)/maxdist)
# dim(pairwise.comps)
# max(pairwise.comps)
net.comps = nrow(pairwise.comps)
nspp = nrow(sampsByCompounds)
ncomps = ncol(sampsByCompounds)
pairwise.spp = as.data.frame(matrix(0,nrow = nspp, ncol = nspp))
names(pairwise.spp) = row.names(sampsByCompounds)
row.names(pairwise.spp) = row.names(sampsByCompounds)
sampsCompsStand = sampsByCompounds

### make this into a binary matrix of metabolite dissimilarity (presence/absence matrix)
sampsCompsStand=as.data.frame(as.matrix((sampsCompsStand>0)+0))


for(i in 1:nrow(sampsByCompounds)){	
  sampsCompsStand[i,] = sampsByCompounds[i,]/sum(sampsByCompounds[i,])
}

diags = pairwise.spp
for (k in 1:nspp){
  sppX = as.character(row.names(sampsCompsStand)[k])
  cat("Comparing ", sppX, " to itself", "\n", sep = "")
  sppXonly = sampsCompsStand[k,which(sampsCompsStand[k,]>0)]
  ncomps = length(sppXonly)
  pairwise.comps.temp = pairwise.comps[names(sppXonly),names(sppXonly)]
  diags[k,k] = sum(((outer(as.numeric(sppXonly), as.numeric(sppXonly)))*pairwise.comps.temp),na.rm = T)	
}
# save(sampsByCompounds, pairwise.comps, diags, file = outfile)
for (i in 1:nspp){
  spp1 = as.character(row.names(sampsCompsStand)[i])
  for (j in i:nspp){
    spp2 = as.character(row.names(sampsCompsStand)[j])
    cat("Comparing ", spp1, " to ", spp2, "\n", sep = "")
    #identify which compounds are in each species
    spp1comps = sampsCompsStand[spp1,]
    spp2comps = sampsCompsStand[spp2,]
    spp_pair = rbind(spp1comps,spp2comps)
    paircomps = spp_pair[,which(colSums(spp_pair)>0)]
    #make a pairwise.comps matrix for only those compounds found in either species in the species pair
    ncomps = ncol(paircomps)
    pairwise.comps.temp = pairwise.comps[names(paircomps),names(paircomps)]
    pairwise.spp[i,j] = pairwise.spp[j,i] = sum(((outer(as.numeric(paircomps[1,]), as.numeric(paircomps[2,])))*pairwise.comps.temp), na.rm = T)/max(diags[i,i], diags[j,j])
  }
}
cscs = pairwise.spp
# save(sppByCompounds, sampsByCompounds, pairwise.comps, sampsCompsStand, diags, cscs, heat.def, siri.itol, tree, meta, feat, file = outfile)

cat("Completed cacluation of CSCS for all sample pairs","\n")

### Bray-Curtis distance matrix of metabolites 
braydist = as.data.frame(as.matrix(vegdist(sampsByCompounds, diag = T, upper = T)))

row.names(braydist) = gsub(".mzXML.Peak.area", "", row.names(braydist))
row.names(braydist) = gsub("X", "", row.names(braydist))
names(braydist) = gsub(".mzXML.Peak.area", "", names(braydist))
names(braydist) = gsub("X", "", names(braydist))

### BINARY (presence/absence) bray-curtis distance matrix of metabolites (diagonals = 0))
braydist.bin = as.data.frame(as.matrix(vegdist(sampsByCompounds, diag = T, upper = T,binary=T)))

row.names(braydist.bin) = gsub(".mzXML.Peak.area", "", row.names(braydist.bin))
row.names(braydist.bin) = gsub("X", "", row.names(braydist.bin))
names(braydist.bin) = gsub(".mzXML.Peak.area", "", names(braydist.bin))
names(braydist.bin) = gsub("X", "", names(braydist.bin))

braydist[1:5,1:5]
braydist.bin[1:5,1:5]

comm.mat = as.data.frame(matrix(1, nrow=2, ncol = nrow(heat.class)))
names(comm.mat) = row.names(heat.class)
qtree.class = prune.sample(samp = comm.mat, phylo = qtree.itol)

heat.class = heat.class[qtree.class$tip.label,]

#### dataframe of Faith's metabolomic distance and metabolite richness for each samples
pd.class = pd(t(heat.class), tree = qtree.class)
row.names(pd.class) = gsub(".mzXML.Peak.area", "", row.names(pd.class))
row.names(pd.class) = gsub("X", "", row.names(pd.class))


save(braydist, braydist.bin, pd.class, class, date, meta,
     file = outfile)
# cscs



#### After CSCS and Bray-Curtis, show heat map, hierarchical clustering, and NMDS plot, then maybe LASSO

#### remove species pools from heat.real
#### fix CSCS calculation


########################
########################
######################## Heat map
########################
########################

heatmap.ht = heat.class
for(i in 1:ncol(heatmap.ht)){
  names(heatmap.ht)[i] = strsplit(names(heatmap.ht)[i], "[.]")[[1]][1]
}

heatmap = heatmap(as.matrix(log10(heatmap.ht + 1)))


########################
########################
######################## NMDS plot, show that conspecific individuals are most chemically similar
########################
########################



library(MASS)
library(viridis)
colors = viridis(7)

nmds.samps = isoMDS(as.dist(braydist), k=2)
nmds.samps.points = nmds.samps$points

ColorCode = rep("black", ncol(braydist))
ColorCode[grep("CARYCO", names(braydist))] = colors[2]
ColorCode[grep("CARYGL", names(braydist))] = colors[3]
ColorCode[grep("CARYO2", names(braydist))] = colors[4]
ColorCode[grep("CARYTE", names(braydist))] = colors[5]
ColorCode[grep("CARYTO", names(braydist))] = colors[6]
ColorCode[grep("JUGLNI", names(braydist))] = colors[7]

plot(nmds.samps.points[,1],nmds.samps.points[,2], 
     col = ColorCode, pch = 16, cex = 2, xlab = "NMDS Axis 1", ylab = "NMDS Axis 2", 
     cex.axis = 1.5, cex.lab = 1.5)


########################
########################
######################## "metabolome-wide association study" using LASSO regression in 'glmnet'
########################
########################




########################
########################
######################## Differential expression analysis using 'deseq2'
########################
########################

