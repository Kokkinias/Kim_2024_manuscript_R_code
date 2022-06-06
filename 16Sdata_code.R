#demonstration R script to manage data and to create 16S rRNA figures/analysis

library(readxl) #to read excel files 
library(ggplot2) #to make plots 
library(vegan) # for ecology tests
library(dplyr) # for data management
library(tidyr) # for data management 
library(stats)# for statistics 

setwd("/pathtodatafiles/")

counts <- read.csv("16S_featuretable.csv",header = TRUE, row.names = "ID")
counts=as.data.frame(counts)

#Filter
#ASV should have at least 10 reads and be in at least 5 samples 
filtered.counts <- counts[rowSums(counts>10)> 5, ] 
write.csv(as.data.frame(filtered.counts), 
          file="16S_filtered_10r5s_featuretable.csv")

#Calculate Relative Abundance
data <- lapply(filtered.counts[], function(x) x/sum(x))
data= as.data.frame(data)
write.csv(as.data.frame(data), 
          file="16S_filtered_10r5s_relabun.csv")

#Bringing in feature table of rel abun
#Excel file includes metadata
#features as columns, samples as rows
features = read_excel("16S_filtered_10r5s_relabun.xlsx", sheet="NMDS_10r5s", col_names = TRUE)
colnames(features) <- c ( "ID","day","highorlow","both")

#make a metadata frame
df <- features [,5:ncol(features)]
df = as.data.frame(sapply(df, as.numeric))
log_df = log(df+1)

########
#NMDS 
########
#LOG transformation
set.seed(3)
bray_dist = metaMDSdist(log_df, distance = "bray", 
                        autotransform = FALSE, noshare = 0.2, trace = 1)

NMDS_Bray <-metaMDS(log_df, distance = "bray",
                    autotransform = FALSE, maxit=800, noshare = 0.2, trace = 1)
#log stress= 0.165

data.scrs = as.data.frame(scores(NMDS_Bray))

#add columns to data frame 
data.scrs$highorlow = features$highorlow
data.scrs$day = features$day
data.scrs$both = features$both
data.scrs$ID = features$ID
data.scrs$day <- as.factor(data.scrs$day)

## Create an NMDS plot with ggplot
plot_getmm = ggplot(data.scrs, aes(x=NMDS1, y=NMDS2)) +
  scale_color_manual(values=c('#da3832','#2f318d','#5bc0eb','#76b956'))+
  geom_point(size= 3, aes(color=day))
plot_getmm

#NMDS Statistics 
##########Day SIGNIFICANT
##anosim and mrpp
anosim <- anosim(bray_dist, features$day, permutations = 999)
anosim
#ANOSIM statistic R: 0.2965 
#Significance: 0.001 

mrpp_data = mrpp(bray_dist, features$day)
mrpp_data
#Chance corrected within-group agreement A: 0.1013 
#Based on observed delta 0.5536 and expected delta 0.6161 
#Significance of delta: 0.001 

#high or low salmonella response 
##########highorlow SIGNIFICANT
##anosim and mrpp
anosim <- anosim(bray_dist, features$highorlow, permutations = 999)
anosim
#ANOSIM statistic R: 0.2876 
#Significance: 0.001

mrpp_data = mrpp(bray_dist, features$highorlow)
mrpp_data
#Chance corrected within-group agreement A: 0.03185 
#Based on observed delta 0.5964 and expected delta 0.6161 
#Significance of delta: 0.003

#Beta dispersion
#Import rel abun table without Salmonella 
features = read_excel("16S_filtered_10r5s_relabun.xlsx", sheet="NMDS_10r5s_nosal", col_names = TRUE)
colnames(features) <- c ( "ID","day","highorlow","both")
#make a metadata frame
df <- features [,5:ncol(features)]
df = as.data.frame(sapply(df, as.numeric))

#Euclidean  
dst <- dist(df)
data.bd <- betadisper(dst,features$both)
data.bd
#Stats for beta dispersion
anova(data.bd)
permutest(data.bd, pairwise = TRUE)

######PROCRUSTES
#MICROBES
#make a metadata frame for microbes
features = read_excel("16S_filtered_10r5s_relabun.xlsx", sheet="NMDS_10r5s", col_names = TRUE)
colnames(features) <- c ( "ID","day","highorlow","both")

df <- features [,5:ncol(features)]
df = as.data.frame(sapply(df, as.numeric))
#PCA for 16S microbes 
pca.micro<-prcomp(df, center=TRUE, scale.=TRUE)

#LIPIDS
lipids_pos = read_excel("16S_filtered_10r5s_relabun.xlsx", sheet="lipid_positive", col_names = TRUE)
colnames(lipids_pos) <- c ( "ID","day","highorlow","both")

lipids_neg = read_excel("16S_filtered_10r5s_relabun.xlsx", sheet="lipid_negative", col_names = TRUE)
colnames(lipids_neg) <- c ( "ID","day","highorlow","both")

#make a metadata frame for lipids in positive mode
df <- lipids [,5:ncol(lipids_pos)]
df = as.data.frame(sapply(df, as.numeric))
#remove all zero columns
nozerocol <- df[, colSums(df != 0) > 0]
#PCA for lipids in positive mode
pca.lipid1<-prcomp(nozerocol, center=TRUE, scale.=TRUE)

#make a metadata frame for lipids in negative mode
df2 <- lipids [,5:ncol(lipids_neg)]
df2 = as.data.frame(sapply(df2, as.numeric))
#remove all zero columns
nozerocol2 <- df2[, colSums(df2 != 0) > 0]
#PCA for lipids in positive mode
pca.lipid2<-prcomp(nozerocol2, center=TRUE, scale.=TRUE)

#PROCRUSTES
pro <- procrustes(X = pca.micro, Y = pca.lipid1, symmetric = TRUE)
pro
#Positive Procrustes sum of squares:0.5531
pro <- procrustes(X = pca.micro, Y = pca.lipid2, symmetric = TRUE)
pro
#Neg Procrustes sum of squares:0.5507

protest(X = pca.micro, Y = pca.lipid1, scores = "sites", permutations = 999)
#Positive 
#Procrustes Sum of Squares (m12 squared):        0.5531 
#Correlation in a symmetric Procrustes rotation: 0.6685 
#Significance:  0.288

protest(X = pca.micro, Y = pca.lipid2, scores = "sites", permutations = 999)
#Negative
#Procrustes Sum of Squares (m12 squared):        0.5507 
#Correlation in a symmetric Procrustes rotation: 0.6703
#Significance:  0.149
