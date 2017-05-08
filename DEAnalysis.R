# DIFFERENTIAL EXPRESSION ANALYSIS USING LIMMA

#source("http://bioconductor.org/biocLite.R")
#biocLite("limma")
library("limma")
library("ggplot2")
library("ggfortify")
library("RColorBrewer")

setwd("/home/arubio/Documents/PROJECTS/040417-Angel_Proteomics/")

# READ DATA -------------------------------------------------------------
# Creating the Expression Set
# Read expression file as Matrix
exprs <- as.matrix(read.table("./RAW/data_zq_subset_renamed.txt", header=TRUE, sep = "\t",
                              row.names = 1,
                              as.is=TRUE))
# For some reason there is an extra column
exprs <- exprs[,1:8]
colnames(exprs) <- c("WT.a", "WT.b","D14.a","D14.b", "D14.c", "D21.a", "D21.b", "D21.c")
head(exprs)
#Features  Samples 
#2657        8

# PCA of the data
# features in cols, observations in rows
t <- t(exprs)
# Add info about the groups
groups <- data.frame(c(rep("WT", 2), rep("D14", 3), rep("D21", 3)))
colnames(groups) <- c("group")
bind <- cbind(t,groups)
# Calculate PCA
pca <- prcomp(t)
# Plot
autoplot(pca, label= TRUE, data=bind, colour='group')
# Save to pdf
pdf(file="./RESULTS/PCA.pdf")
autoplot(pca, label= TRUE, data=bind, colour='group')
dev.off()

# read phenotipic data file
pData <- read.table("./ANALYSIS/01.Differential_Expression/All_zq_subset_phenotipic_data.txt",row.names=1, header=TRUE, sep="\t")
#pData
# TYPE TIME
# WT_127_N        WT  D21
# WT_127_C        WT  D21
# D14_128_N CD196CTR  D14
# D14_128_C CD196CTR  D14
# D14_129_N CD196CTR  D14
# D21_129_C CD196CTR  D21
# D21_130_N CD196CTR  D21
# D21_130_C CD196CTR  D21

# Check if pheno data and expression data info matches
all(rownames(pData)==colnames(exprs))
# TRUE

# Now we can add metadata to the pData data frame
metadata <- data.frame(labeldescription=c("Mouse type", "Day after experiment"))
# Costruct the Annotated Data Frame using the Pdata data frame and the metadata
phenoData <- new ("AnnotatedDataFrame", data=pData, varMetadata=metadata)
# Create the expression set with it's phenotype data
eset <- ExpressionSet(assayData=exprs, phenoData=phenoData)

# READ TARGETS -------------------------------------------------------------
# What has been tested in each color of the array; so the A and B in log2(A/B).
# Name of the channels must be Cy3 and Cy5 for limma to be able to create model matrix.
# In this case we tested each protein (A) with a technical reference (B).
targets <- readTargets("./ANALYSIS/01.Differential_Expression/targets.txt", sep="")
# SampleID Cy3 Cy5
# 1  WT_127_N Ref  WT
# 2  WT_127_C Ref  WT
# 3 D14_128_N Ref D14
# 4 D14_128_C Ref D14
# 5 D14_129_N Ref D14
# 6 D21_129_C Ref D21
# 7 D21_130_N Ref D21
# 8 D21_130_C Ref D21

# CREATE MODEL MATRIX
design <- modelMatrix(targets,ref="Ref")
# Fitting linear model
# Calculate statistics with the targets info. (Replicates, reference, etc)
fit <- lmFit(eset, design)

# CREATE CONTRAST MATRIX
# Determines which comparisons we want to make. We want to test each condition to each other:
# So D14-WT=log2(D14/WT), etc
cont.matrix <- makeContrasts("D14/WT"="D14-WT","D21/WT"="D21-WT","D21/D14"="D21-D14",levels=design)

# E-BAYES -------------------------------------------------------------------
# Statistics method adjusted for a small N and that accounts for global variance. 
# (ANOVA wouldnt be possible with this N size)
fit2  <- contrasts.fit(fit, cont.matrix)
fit2  <- eBayes(fit2)

# GET D.E. PROTEINS ---------------------------------------------------------
# Get a knack of the results
# colnames(fit2)
# topTable(fit2,coef=1,adjust="fdr")
# topTable(fit2,coef=2,adjust="fdr")
# topTable(fit2,coef=3,adjust="fdr")

# The outcome of each hypothesis test can be assigned using
results <- decideTests(fit2)

# Visualize data
vennDiagram(results)

# Getting more info from the diagram:
# Upregulated
vennCounts(results,include="up")
# D14/WT D21/WT D21/D14 Counts
# 1      0      0       0   2411
# 2      0      0       1     60
# 3      0      1       0     63
# 4      0      1       1    100
# 5      1      0       0      0
# 6      1      0       1      0
# 7      1      1       0     13
# 8      1      1       1     10
vennDiagram(vennCounts(results,include="up"))
title("Upregulated", line = -1)

# Downregulated
vennCounts(results,include="down")
# D14/WT D21/WT D21/D14 Counts
# 1      0      0       0   2492
# 2      0      0       1     45
# 3      0      1       0     51
# 4      0      1       1     58
# 5      1      0       0      7
# 6      1      0       1      0
# 7      1      1       0      4
# 8      1      1       1      0
vennDiagram(vennCounts(results,include="down"))
title("Downregulated", line = -1)


#"D14/WT"=1,"D21/WT"="D21-WT"2,"D21/D14"=3
# Get proteins different in D21/WT but not in the other contrasts
D21_WT_prots<-which(results[,2]!=0 & results[,1]==0 & results[,3]==0)
write.table(results[D21_WT_prots,], file="./RESULTS/D21_WT_prots.xls",sep="\t",quote=FALSE, col.names=NA)
# Get proteins different in D14/WT but not in the other contrasts
D14_WT_prots<-which(results[,1]!=0 & results[,2]==0 & results[,3]==0)
write.table(results[D14_WT_prots,], file="./RESULTS/D14_WT_prots.xls",sep="\t",quote=FALSE, col.names=NA)
# Get proteins different in D21/14 but not in the other contrasts
D21_14_prots<-which(results[,3]!=0 & results[,2]==0 & results[,1]==0)
write.table(results[D21_14_prots,], file="./RESULTS/D21_14_prots.xls",sep="\t",quote=FALSE, col.names=NA)

# Get proteins different in D21/14 and D14/WT but not in D21/WT
D21_14_D14_WT_prots<-which(results[,3]!=0 & results[,2]==0 & results[,1]!=0)
write.table(results[D21_14_D14_WT_prots,], file="./RESULTS/D21_14_D14_WT_prots.xls",sep="\t",quote=FALSE, col.names=NA)
# Get proteins different in D21/WT and D14/WT but not in D21/14
D21_WT_D14_WT_prots <- which(results[,1]!=0 & results[,2]!=0 & results[,3]==0)
write.table(results[D21_WT_D14_WT_prots,], file="./RESULTS/D21_WT_D14_WT_prots.xls",sep="\t",quote=FALSE, col.names=NA)
# Get proteins different in D21/14 and D21/WT but not in D14/WT
D21_14_D21_WT_prots <- which(results[,1]==0 & results[,2]!=0 & results[,3]!=0)
write.table(results[D21_14_D21_WT_prots,], file="./RESULTS/D21_14_D21_WT_prots.xls",sep="\t",quote=FALSE, col.names=NA)
# Get proteins different in all three contrasts
all_prots <- which(results[,1]!=0 & results[,2]!=0 & results[,3]!=0)
write.table(results[all_prots,], file="./RESULTS/all_prots.xls",sep="\t",quote=FALSE, col.names=NA)
# 
# # Export all differentially expressed to one file
all_prots <- which(results[,1]!=0 | results[,2]!=0 | results[,3]!=0)
write.table(results[all_prots,], file="./ANALYSIS/01.Differential_Expression/fit_DE.txt",sep="\t",quote=FALSE, col.names=NA)

# 
# # Or write one file per contrast
# numProts <- nrow(exprs)
# completeTableD14_WT <- topTable(fit2,coef=1,number=numProts, adjust="fdr")
# write.table(completeTableD14_WT,file="./RESULTS/D14_WT.xls",sep="\t",quote=FALSE, col.names=NA)
# completeTableD21_WT <- topTable(fit2,coef=2,number=numProts, adjust="fdr")
# write.table(completeTableD21_WT,file="./RESULTS/D21_WT.xls",sep="\t",quote=FALSE, col.names=NA)
# completeTableD21_14 <- topTable(fit2,coef=3,number=numProts, adjust="fdr")
# write.table(completeTableD21_14,file="./RESULTS/D21_D14.xls",sep="\t",quote=FALSE, col.names=NA)

# HEATMAPS
# color palette
hmcol<-brewer.pal(11,"RdBu")
# Heatmap of all the dataset

heatmap(exprs, col=hmcol)
# Heatmap of all the 405 differentially expressed genes
de.results <- which(results[,1]!=0 | results[,2]!=0 | results[,3]!=0)
de.prots <- results[de.results,]
myprots <- rownames(de.prots)
# Select subset from original file
subset <- exprs[myprots,]
# Create the actual heatmap
heatmap(subset, col=hmcol)
# save to pdf
pdf(file="./RESULTS/heatmap_all_de_prots.pdf")
heatmap(subset, col=hmcol)
dev.off()

