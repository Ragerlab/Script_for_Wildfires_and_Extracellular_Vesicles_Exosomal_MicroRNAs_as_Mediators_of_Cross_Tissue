#################################################################################################
#################################################################################################
#### Analysis of RNA-Seq data for the Wildfire Project
####
#### Using DESeq2 for sample normalization and statistical analysis
####
#### Because this is just a toxicity analysis, code does not include covariate imputation, and just an optional SVA step
#### 
#### 
#### Code drafted by Alexis Payton, Julia Rager, Lauren Koval, & Vennela Avula
#### Lasted updated: February 28, 2021
#################################################################################################
#################################################################################################


sessionInfo()
rm(list=ls())

#################################################################################################
#################################################################################################
#### Installing appropriate packages in R (if already have these installed, SKIP THIS STEP)
#################################################################################################
#################################################################################################
# install.packages("data.table")

# Note that BiocManager allows you to access the Bioconductor repo in which DeSeq2 package is. DeSeq2 pacakge has a "tibble" "RcppArmadillo" and "rlang" dependency.
# to install RcppArmadillo, you must install a Fortan compiler (I used gfortran-6.1.pkg, https://cran.r-project.org/bin/macosx/tools/).

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# install.packages("tibble")
# install.packages("rlang")
# install.packages("RcppArmadillo")
# install.packages("janitor")

# Install DESeq and SVA: 
# BiocManager::install("DESeq2", version = "3.11")
# BiocManager::install("sva", version ="3.11")
# BiocManager::install("BiocManager")
# BiocManager::install("edgeR")

# Install RUVseq
# BiocManager::install("RUVSeq")

# Install plot packages
# install.packages("ggbeeswarm") # note that I had to manually install this package using tools --> install packages
# install.packages("gridExtra")
# install.packages("tidyverse") # Install tidyverse packages 
# install.packages("gplots")
# install.packages("RColorBrewer")
# install.packages("scales") # Install scales (needed for tidyverse activation)
# install.packages("plotly")




#################################################################################################
#################################################################################################
#### Activating appropriate packages in R
#################################################################################################
#################################################################################################

# Activating appropriate libraries
library(data.table)
library(rlang)
library(DESeq2)
library(limma)
library(sva)
library(ggplot2)
library(stats)
library(scales)
library(ggbeeswarm)
library(gridExtra)
library(tidyverse)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(edgeR) 
library(plotly)
library(RUVSeq)
library(janitor)

#################################################################################################
#################################################################################################
#### Set working directory, output folder path, etc
#################################################################################################
#################################################################################################

# Set working directory
setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Wildfire/Experiments/Experiment_2/1_Input")
getwd()

# Create an output folder (make sure to make the folder first, and then point to it here)
Output <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Experiment_2/2_Output"
Output_StatResults <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Experiment_2/2_Output/StatResults"
Output_Figures <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Experiment_2/3_Figures"
Output_RUV <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Experiment_2/2_Output/StatResults_RUV"
cur_date <- "022821"

#################################################################################################
#################################################################################################
#### Loading count data and metadata (sample information data)
#### Organizing and filtering to keep samples that are present all files
#################################################################################################
#################################################################################################


# Read in count data
countdata <- read.csv(file = 'SP0222_gene_counts.csv', check.names = FALSE)
dim(countdata)
# 171 samples, 30146 mRNAs

# visualize this data quickly by viewing top left corner:
countdata[1:3,1:6]

# Check for duplicate mRNA IDs in the countdata 
Dups <- duplicated(countdata[,1])
summary(Dups)
# Not an issue for this wildfire dataset


# Read in metadata (sample information file)
subjectinfo <- read.csv(file = "Sample_Info_112520.csv", check.names = FALSE)
dim(subjectinfo) #170 rows, 9 col

# Visualize this data quickly by viewing top left corner:
subjectinfo[1:3,1:7]

#reading in table of contrasts
contrasts = read_csv("Table_of_Contrasts_121120.csv")
head(contrasts)

#################################################################################################
#################################################################################################
#### Create dataframes that are formatted for proceeding code, as well as DESeq2 functions
#### countdata and coldata
#################################################################################################
#################################################################################################

# Creating the coldata object, based on information in the subjectinfo file
coldata <- subjectinfo

# Set the rownames of coldata and column names of countdata to be in the same order 
countdata <- setcolorder(countdata, as.character(coldata$SampleID_BioSpyderCountFile))

# Make Gene the row name so it isn't renamed with the sample ID
countdata <- countdata %>% column_to_rownames("Gene")

# Replace the countdata column names (biospyder ids) with the treatment ids
colnames(countdata) <- coldata$ID

# Double checking that the same variables appear between the two dataframes
setequal(as.character(coldata$ID), colnames(countdata))

# Additionally checking that not only the sets of variables are the same, but that they are in the same order
identical(as.character(coldata$ID), colnames(countdata))


###################################################################################################
###################################################################################################
#### Background Filter
###################################################################################################
###################################################################################################

# Total number of samples
nsamp <- ncol(countdata)

# Median expression level across all genes and all samples
total_median <- median(as.matrix(countdata))

# Add back in the Gene column
countdata <- countdata %>% rownames_to_column("Gene")

# Get list of genes that have an expression greater than the total median in at least 20% of the samples
genes_above_background <- countdata %>% 
  pivot_longer(cols=!Gene, names_to = "sampleID", values_to="expression") %>% 
  mutate(above_median=ifelse(expression>total_median,1,0)) %>% 
  group_by(Gene) %>% 
  summarize(total_above_median=sum(above_median)) %>% 
  filter(total_above_median>.2*nsamp) %>% 
  select(Gene)

# filter countdata for only the genes above background
countdata <- left_join(genes_above_background, countdata, by="Gene")


# Write out filtered countdata
write.csv(countdata, paste0(Output,"/", cur_date, "_RawCounts_AboveBack.csv"), row.names=TRUE)

#################################################################################################
#################################################################################################
#### Subject Filter
#################################################################################################
#################################################################################################

# Transpose filtered countdata and confirm all samples have expression
countdata_T <- countdata %>% 
  pivot_longer(cols=!Gene, names_to="sampleID",values_to="expression") %>% 
  pivot_wider(names_from=Gene, values_from=expression) %>% 
  adorn_totals(where="col", name="rowsum")
 
# all samples have some expression
nrow(countdata_T %>% filter(rowsum==0))



#################################################################################################
#################################################################################################
#### RNASeq QA/QC on raw count data to identify potential outlier samples
#### This may or may not result in another filter step
#################################################################################################
#################################################################################################
# One way to evaluate sequencing data quality / identify potential outliers is through Principal Component Analysis (PCA)
# PCA helps in identifying outlying samples for quality control, and gives a feeling for the principal causes of variation in a dataset

countdata <- countdata %>% column_to_rownames("Gene")

# Calculate principal components using transposed count data
pca <- prcomp(t(countdata))
screeplot(pca)
# this scree plot indicates that nearly all variation is explained in PC1,2 and 3, 
# thus PCA is a useful exercise for this dataset 


# Make dataframe for PCA plot generation using first three components
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], Sample=colnames(countdata))
# Add attributes (covariates from coldata) to pca df
pca_df <- merge(pca_df, coldata, by.x="Sample", by.y="ID")


# Calculating percent of the variation that is captured by each principal component
pca_percent <- round(100*pca$sdev^2/sum(pca$sdev^2),1)


# Generating PCA plot annotated by plate batch 

# 2D
ggplot(pca_df, aes(PC1,PC2, color = PlateBatch))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output_Figures,"/", cur_date,"_PCA_PlateBatch.png"))
dev.off()


# Generating PCA plot annotated by treatment 

#2D
ggplot(pca_df, aes(PC1,PC2, color = Treatment))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output_Figures,"/", cur_date,"_PCA_Treatment.png"))
dev.off()


# Generating PCA plot annotated by tissue

# 2D
ggplot(pca_df, aes(PC1,PC2, color = Tissue))+
  geom_point(size=6) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output_Figures,"/", cur_date,"_PCA_Tissue.png"))
dev.off()


# all the heart samples seem to be grouped together, whcih makes sense
# lung samples display more variance than heart - so we may want to capture unwanted variance through SVA or RUV techniques


# Lets identify which samples are which

# 2D
ggplot(pca_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output_Figures,"/", cur_date,"_PCA_IDsLabeled.png"))
dev.off()


# This would suggest that M112_RedOakFlame may be a potential outlier, but we can investigate further
# We will conduct heirachical clustering to investigate further and to determine if we should remove these samples


# create a dataframe with samples as rows and genes as columns to input into hclust
countdata_forclustering <- t(countdata)
countdata_forclustering[1:5,1:10]

# run hierarchical clustering
dev.off()
hc <-hclust(dist(countdata_forclustering))
jpeg(filename= paste0(Output_Figures,"/",cur_date, "_Hierachical_clustering_outliers.jpeg"),width = 4500, height = 350)
plot(hc)
dev.off()


# Here, the one aforementioned sample (M112_RedOakFlame) was separate from the rest, so should be removed



#################################################################################################
#################################################################################################
#### Removing sample(s) with potential QA/QC issues
#################################################################################################
#################################################################################################


# First, pull all the IDs from the count data
keep_samples <- colnames(countdata)
# Note that in this case, this represents all 170 IDs

# Then, remove the one sample outlier from the "keep_samples" vectors
keep_samples <- keep_samples[! keep_samples %in% c("M112_RedOakFlame")]
# Note that this represents all remaining 169 IDs

# Filter the countdata dataframe for these IDs
countdata <- countdata[, colnames(countdata) %in% keep_samples]

# Filter the coldata dataframe for these IDs
coldata <- coldata[coldata$ID %in% keep_samples, ]


# Export the raw data for just the included samples to potentially generate plots outside of R
write.csv(countdata, paste0(Output,"/", cur_date, "_RawCounts_SamplesIncluded.csv"), row.names= TRUE)
write.csv(coldata, paste0(Output,"/", cur_date, "_SampleInfo_SamplesIncluded.csv"), row.names= FALSE)





#################################################################################################
#################################################################################################
#### Setting up the actual DESeq2 function and associated algorithm ("experiment")
#################################################################################################
#################################################################################################
# More information can be found about the DeSeq2 package at https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8, as well as many other resources, as this package is well documented


# Specify contrasts

#design to eliminate the intercept
design <- model.matrix(~0 + Group, data=coldata)

#set up comparisons
#whatever is first in each statement is the numerator in the fold change estimate (treatment - control is like saying treatment/control)
#putting a plus sign in the formula looks at the additive affect of strain on the exposure effect; minus sign see the interaction of strain and exposure:
my_contrasts <- makeContrasts(RedOakFlame_24h_Heart = GroupRedOakFlame_24h_Heart - GroupSaline_24h_Heart,
                              PeatFlame_24h_Heart = GroupPeatFlame_24h_Heart - GroupSaline_24h_Heart, 
                              RedOakSmolder_24h_Heart = GroupRedOakSmolder_24h_Heart - GroupSaline_24h_Heart, 
                              PeatSmolder_24h_Heart = GroupPeatSmolder_24h_Heart - GroupSaline_24h_Heart, 
                              RedOakSmolder_4h_Lung = GroupRedOakSmolder_4h_Lung - GroupSaline_4h_Lung,
                              PeatSmolder_4h_Lung = GroupPeatSmolder_4h_Lung - GroupSaline_4h_Lung,
                              PineNeedlesSmolder_4h_Lung = GroupPineNeedlesSmolder_4h_Lung - GroupSaline_4h_Lung,
                              PineSmolder_4h_Lung = GroupPineSmolder_4h_Lung - GroupSaline_4h_Lung,
                              EucalyptusSmolder_4h_Lung = GroupEucalyptusSmolder_4h_Lung - GroupSaline_4h_Lung,
                              RedOakFlame_4h_Lung = GroupRedOakFlame_4h_Lung - GroupSaline_4h_Lung,
                              PeatFlame_4h_Lung = GroupPeatFlame_4h_Lung - GroupSaline_4h_Lung,
                              PineNeedlesFlame_4h_Lung = GroupPineNeedlesFlame_4h_Lung - GroupSaline_4h_Lung,
                              PineFlame_4h_Lung = GroupPineFlame_4h_Lung - GroupSaline_4h_Lung,
                              EucalyptusFlame_4h_Lung = GroupEucalyptusFlame_4h_Lung - GroupSaline_4h_Lung,
                              LPS_4h_Lung = GroupLPS_4h_Lung - GroupSaline_4h_Lung,
                              RedOakSmolder_24h_Lung = GroupRedOakSmolder_24h_Lung - GroupSaline_24h_Lung,
                              PeatSmolder_24h_Lung = GroupPeatSmolder_24h_Lung - GroupSaline_24h_Lung,
                              PineNeedlesSmolder_24h_Lung = GroupPineNeedlesSmolder_24h_Lung - GroupSaline_24h_Lung,
                              PineSmolder_24h_Lung = GroupPineSmolder_24h_Lung - GroupSaline_24h_Lung,
                              EucalyptusSmolder_24h_Lung = GroupEucalyptusSmolder_24h_Lung - GroupSaline_24h_Lung,
                              RedOakFlame_24h_Lung = GroupRedOakFlame_24h_Lung - GroupSaline_24h_Lung,
                              PeatFlame_24h_Lung = GroupPeatFlame_24h_Lung - GroupSaline_24h_Lung,
                              PineNeedlesFlame_24h_Lung = GroupPineNeedlesFlame_24h_Lung - GroupSaline_24h_Lung,
                              PineFlame_24h_Lung = GroupPineFlame_24h_Lung - GroupSaline_24h_Lung,
                              EucalyptusFlame_24h_Lung = GroupEucalyptusFlame_24h_Lung - GroupSaline_24h_Lung,
                              LPS_24h_Lung = GroupLPS_24h_Lung - GroupSaline_24h_Lung,
                              levels=design)



# Ensuring that the appropriate variable types are recognized within the DESeq2 algorithm (e.g., factor vs numeric)
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Tissue <- as.factor(coldata$Tissue)
coldata$Group <- as.factor(coldata$Group)




#########################################################
# Set up DESeq2 experiment
#########################################################
dds = DESeqDataSetFromMatrix(countData = countdata,
                             colData = coldata,
                             design = ~0+Group) #design matches the design above for specifying contrasts


dds <- estimateSizeFactors(dds)
#sizeFactors(dds) #check size factors

# normalized counts
normcounts <- counts(dds, normalized=TRUE)
write.csv(normcounts, paste0(Output,"/",cur_date, "_NormCounts.csv"), row.names=TRUE)


# pseudocounts
ps_normcounts <- normcounts + 1
write.csv(ps_normcounts, paste0(Output,"/",cur_date, "_NormCounts_ps.csv"),row.names=TRUE)


# log2 pseudocounts (y=log2(n+1))
log2normcounts <- log2(normcounts+1)
write.csv(log2normcounts, paste0(Output,"/",cur_date, "_NormCounts_pslog2.csv"), row.names=TRUE)


dds <- DESeq(dds, betaPrior=FALSE)


# loop through and write csvs of results for all contrasts
count_list <- c()
n <-1 
for (contrast in colnames(my_contrasts)){
  cat(contrast)
  cat("\n")
  res <- results(dds, pAdjustMethod = "BH", contrast = my_contrasts[,contrast])  #Statistical output with multiple test correction by the default, BH (aka FDR)
  summary(res)
  ordered <- as.data.frame(res[order(res$padj),])
  gene_count_padj <- dim(ordered %>% filter(padj<0.1))[1]
  write.csv(ordered, paste0(Output_StatResults,"/", cur_date, "_Stat_Results_", contrast ,".csv"))
  count_list[n] <- gene_count_padj 
  n <- n+1
  cat("\n\n")
}

# table with the number of significantly differentially expressed Genes for each contrast at padj= 0.1
counts_df <- data.frame(Group=colnames(my_contrasts),
                        Count=count_list)
write.csv(counts_df, paste0(Output,"/", cur_date, "_Count_Sig_Genes_Table.csv"), row.names = FALSE)

# Variance stabilizing matrix
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd_matrix <-as.matrix(assay(vsd))
write.csv(vsd_matrix, paste0(Output,"/", cur_date, "_VSDCounts.csv"), row.names=TRUE)





#################################################################################################
#################################################################################################
#### RUVs normalization - accounts for variation that occurs within replicate groups
#################################################################################################
#################################################################################################

# following these resources: http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf , https://www.nature.com/articles/nbt.2931 

genes <- rownames(countdata)
x <- as.factor(coldata$Group) #creating a vector of the Groups 

set <- newSeqExpressionSet(as.matrix(countdata),
                           phenoData = data.frame(x, row.names=colnames(countdata))) #storing the data in an object of S4 class SeqExpressionSet to make use of the plotting and normalization functionality of RUVSeq 

# exploratory plots 
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

#construct a matrix specifying the replicates
differences <- makeGroups(x)

set3 <- RUVs(set, genes, k=3, differences) #using k=3 factors of unwanted variation as 2-3 factors are generally enough (reference:David Risso, https://www.nature.com/articles/nbt.2931)
# This returns two pieces of information
##### 1. estimated factors of unwanted variation (phenoData in set)
##### 2. normalized counts obtained by regressing the original counts on the unwanted factors (normalizedCounts in set)



#duplicate coldata
coldata_RUV <-data.frame(coldata)

#add in surrogate W_1
coldata_RUV$group_check <- set3$x
coldata_RUV$W_1 <- set3$W_1


#make sure groups are aligned appropriately
identical(coldata_RUV$Group,coldata_RUV$group_check)


#new design including factors of unwanted variation
design_RUV <- model.matrix(~0 + Group+W_1, data=coldata_RUV)

#set up comparisons
#whatever is first in each statement is the numerator in the fold change estimate (treatment - control is like saying treatment/control)
#putting a plus sign in the formula looks at the additive affect of strain on the exposure effect; minus sign see the interaction of strain and exposure:
my_contrasts_RUV <- makeContrasts(RedOakFlame_24h_Heart = GroupRedOakFlame_24h_Heart - GroupSaline_24h_Heart,
                              PeatFlame_24h_Heart = GroupPeatFlame_24h_Heart - GroupSaline_24h_Heart, 
                              RedOakSmolder_24h_Heart = GroupRedOakSmolder_24h_Heart - GroupSaline_24h_Heart, 
                              PeatSmolder_24h_Heart = GroupPeatSmolder_24h_Heart - GroupSaline_24h_Heart, 
                              RedOakSmolder_4h_Lung = GroupRedOakSmolder_4h_Lung - GroupSaline_4h_Lung,
                              PeatSmolder_4h_Lung = GroupPeatSmolder_4h_Lung - GroupSaline_4h_Lung,
                              PineNeedlesSmolder_4h_Lung = GroupPineNeedlesSmolder_4h_Lung - GroupSaline_4h_Lung,
                              PineSmolder_4h_Lung = GroupPineSmolder_4h_Lung - GroupSaline_4h_Lung,
                              EucalyptusSmolder_4h_Lung = GroupEucalyptusSmolder_4h_Lung - GroupSaline_4h_Lung,
                              RedOakFlame_4h_Lung = GroupRedOakFlame_4h_Lung - GroupSaline_4h_Lung,
                              PeatFlame_4h_Lung = GroupPeatFlame_4h_Lung - GroupSaline_4h_Lung,
                              PineNeedlesFlame_4h_Lung = GroupPineNeedlesFlame_4h_Lung - GroupSaline_4h_Lung,
                              PineFlame_4h_Lung = GroupPineFlame_4h_Lung - GroupSaline_4h_Lung,
                              EucalyptusFlame_4h_Lung = GroupEucalyptusFlame_4h_Lung - GroupSaline_4h_Lung,
                              LPS_4h_Lung = GroupLPS_4h_Lung - GroupSaline_4h_Lung,
                              RedOakSmolder_24h_Lung = GroupRedOakSmolder_24h_Lung - GroupSaline_24h_Lung,
                              PeatSmolder_24h_Lung = GroupPeatSmolder_24h_Lung - GroupSaline_24h_Lung,
                              PineNeedlesSmolder_24h_Lung = GroupPineNeedlesSmolder_24h_Lung - GroupSaline_24h_Lung,
                              PineSmolder_24h_Lung = GroupPineSmolder_24h_Lung - GroupSaline_24h_Lung,
                              EucalyptusSmolder_24h_Lung = GroupEucalyptusSmolder_24h_Lung - GroupSaline_24h_Lung,
                              RedOakFlame_24h_Lung = GroupRedOakFlame_24h_Lung - GroupSaline_24h_Lung,
                              PeatFlame_24h_Lung = GroupPeatFlame_24h_Lung - GroupSaline_24h_Lung,
                              PineNeedlesFlame_24h_Lung = GroupPineNeedlesFlame_24h_Lung - GroupSaline_24h_Lung,
                              PineFlame_24h_Lung = GroupPineFlame_24h_Lung - GroupSaline_24h_Lung,
                              EucalyptusFlame_24h_Lung = GroupEucalyptusFlame_24h_Lung - GroupSaline_24h_Lung,
                              LPS_24h_Lung = GroupLPS_24h_Lung - GroupSaline_24h_Lung,
                              levels=design_RUV)


dds_RUV = DESeqDataSetFromMatrix(countData = countdata,
                             colData = coldata_RUV,
                             design = ~0+Group+W_1) #design matches the design above for specifying contrasts


dds_RUV <- estimateSizeFactors(dds_RUV)

# normalized counts
normcounts_RUV <- counts(dds_RUV, normalized=TRUE)
write.csv(normcounts_RUV, paste0(Output,"/",cur_date, "_NormCounts_RUV.csv"), row.names=TRUE)


# pseudocounts
ps_normcounts_RUV <- normcounts_RUV + 1
write.csv(ps_normcounts_RUV, paste0(Output,"/",cur_date, "_NormCounts_ps_RUV.csv"),row.names=TRUE)


# log2 pseudocounts (y=log2(n+1))
log2normcounts_RUV <- log2(normcounts_RUV+1)
write.csv(log2normcounts_RUV, paste0(Output,"/",cur_date, "_NormCounts_pslog2_RUV.csv"), row.names=TRUE)



dds_RUV <- DESeq(dds_RUV, betaPrior=FALSE) 

# loop through and write csvs of results for all contrasts
count_list_RUV <- c()
n <-1 
for (contrast in colnames(my_contrasts_RUV)){
  cat(contrast)
  cat("\n")
  res <- results(dds_RUV, pAdjustMethod = "BH", contrast = my_contrasts_RUV[,contrast])  #Statistical output with multiple test correction by the default, BH (aka FDR)
  summary(res)
  ordered <- as.data.frame(res[order(res$padj),])
  gene_count_padj <- dim(ordered %>% filter(padj<0.1))[1]
  write.csv(ordered, paste0(Output_RUV,"/",cur_date,"_Stat_Results_RUV_", contrast ,".csv"))
  count_list_RUV[n] <- gene_count_padj 
  n <- n+1
  cat("\n\n")
}

# table with the number of significantly differentially expressed Genes for each contrast at padj= 0.1
counts_df_RUV <- data.frame(Group=colnames(my_contrasts_RUV),
                        Count=count_list_RUV)
write.csv(counts_df_RUV, paste0(Output,"/", cur_date, "_Count_Sig_Genes_Table_RUV.csv"), row.names = FALSE)

# Variance stabilizing matrix
vsd_RUV <- varianceStabilizingTransformation(dds_RUV, blind=FALSE)
vsd_matrix_RUV <-as.matrix(assay(vsd_RUV))
write.csv(vsd_matrix_RUV, paste0(Output,"/", cur_date, "_VSDCounts_RUV.csv"), row.names=TRUE)




