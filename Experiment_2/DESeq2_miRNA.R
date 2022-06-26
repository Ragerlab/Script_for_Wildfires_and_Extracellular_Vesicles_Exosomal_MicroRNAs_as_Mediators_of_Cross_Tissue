#################################################################################################
#################################################################################################
#### Analysis of miRNA-Seq data 
####
#### Using DESeq2 for sample normalization and statistical analysis
####
#### 
#### 
#### Code drafted by Vennela Avula 
#### Lasted updated: April 4, 2021
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
setwd("/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Wildfire/Experiments/Experiment_3/Input")
getwd()

# Create an output folder (make sure to make the folder first, and then point to it here)
Output <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Wildfire/Experiments/Experiment_3/Output"
Output_StatResults <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Wildfire/Experiments/Experiment_3/Output/StatResults"
Output_Figures <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Wildfire/Experiments/Experiment_3/Output/Figures"
# Output_RUV <- "/Users/lkoval/IEHS Dropbox/Rager Lab/Lauren_Koval/LK_Lab_Notebook/Projects/Wildfire/Experiments/Experiment_3/Output/RUV"
cur_date <- "040421"

#################################################################################################
#################################################################################################
#### Loading count data and metadata (sample information data)
#### Organizing and filtering to keep samples that are present all files
#################################################################################################
#################################################################################################


# Read in count data from 2 batches
countdata0 <- read.csv(file = 'batch0_mouse_miRNA_count.csv', check.names = FALSE)
countdata1 <- read.csv(file = 'batch1_mouse_miRNA_count.csv', check.names = FALSE)

# Merge both countdata files to include only miRNAs in both
countdata0 <- countdata0 %>% column_to_rownames("miRNA")
countdata1 <- countdata1 %>% column_to_rownames("miRNA")
countdata <- merge(countdata0,countdata1, 0)
countdata <- as.data.frame(countdata)
dim(countdata)
# 30 samples, 237 miRNAs

# visualize this data quickly by viewing top left corner:
countdata[1:3,1:6]

# Check for duplicate miRNA IDs in the countdata 
Dups <- duplicated(countdata[,1])
summary(Dups)
# Not an issue for this wildfire miRNA dataset


# Read in metadata (sample information file)
subjectinfo <- read.csv(file = "microRNAseq_ID_summaries_032821.csv", check.names = FALSE)
dim(subjectinfo) #60 rows, 9 col

# Visualize this data quickly by viewing top left corner:
subjectinfo[1:3,1:7]

#reading in table of contrasts
contrasts = read_csv("Table_of_Contrasts_032821.csv")
head(contrasts)

#################################################################################################
#################################################################################################
#### Create dataframes that are formatted for proceeding code, as well as DESeq2 functions
#### countdata and coldata
#################################################################################################
#################################################################################################

# Creating the coldata object, based on information in the subjectinfo file
coldata <- subjectinfo

# Make miRNA the row name so it isn't renamed with the sample ID
countdata <- countdata %>% column_to_rownames("Row.names")

# Filter colddata to only include sample names from count data
samples <- colnames(countdata)
coldata <- subset(coldata, Label_from_Core %in% samples)

# Set the rownames of coldata and column names of countdata to be in the same order 
countdata <- setcolorder(countdata, as.character(coldata$Label_from_Core))

# Double checking that the same variables appear between the two dataframes
setequal(as.character(coldata$Label_from_Core), colnames(countdata)) #TRUE

# Additionally checking that not only the sets of variables are the same, but that they are in the same order
identical(as.character(coldata$Label_from_Core), colnames(countdata)) #TRUE 


###################################################################################################
###################################################################################################
#### Background Filter
###################################################################################################
###################################################################################################

# Total number of samples
nsamp <- ncol(countdata) #30 

# Median expression level across all genes and all samples
total_median <- median(as.matrix(countdata)) #1

# Add back in the miRNA column
countdata <- countdata %>% rownames_to_column("miRNA")

# Get list of genes that have an expression greater than the total median in at least 20% of the samples
miRNAs_above_background <- countdata %>% 
  pivot_longer(cols=!miRNA, names_to = "sampleID", values_to="expression") %>% 
  mutate(above_median=ifelse(expression>total_median,1,0)) %>% 
  group_by(miRNA) %>% 
  summarize(total_above_median=sum(above_median)) %>% 
  filter(total_above_median>.2*nsamp) %>% 
  select(miRNA) #155  out of 237 

# filter countdata for only the genes above background
countdata <- left_join(miRNAs_above_background, countdata, by="miRNA")


# Write out filtered countdata
write.csv(countdata, paste0(Output,"/", cur_date, "_RawCounts_AboveBack.csv"), row.names=TRUE)

#################################################################################################
#################################################################################################
#### Subject Filter
#################################################################################################
#################################################################################################

# Transpose filtered countdata and confirm all samples have expression
countdata_T <- countdata %>% 
  pivot_longer(cols=!miRNA, names_to="sampleID",values_to="expression") %>% 
  pivot_wider(names_from=miRNA, values_from=expression) %>% 
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

countdata <- countdata %>% column_to_rownames("miRNA")

# Calculate principal components using transposed count data
pca <- prcomp(t(countdata))
screeplot(pca)
# this scree plot indicates that nearly all variation is explained in PC1,2 and 3, 
# thus PCA is a useful exercise for this dataset 


# Make dataframe for PCA plot generation using first three components
pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], Sample=colnames(countdata))
# Add attributes (covariates from coldata) to pca df
pca_df <- merge(pca_df, coldata, by.x="Sample", by.y="Label_from_Core")


# Calculating percent of the variation that is captured by each principal component
pca_percent <- round(100*pca$sdev^2/sum(pca$sdev^2),1)


# Generating PCA plot annotated by plate batch  

# 2D
ggplot(pca_df, aes(PC1,PC2, color = Batch))+
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


# Lets identify which samples are which

# 2D
ggplot(pca_df, aes(PC1,PC2))+
  geom_text(aes(label=Sample, size=3)) +
  labs(x=paste0("PC1 (",pca_percent[1],"%)"), y=paste0("PC2 (",pca_percent[2],"%)"))
dev.copy(png,paste0(Output_Figures,"/", cur_date,"_PCA_IDsLabeled.png"))
dev.off()

# We will conduct heirachical clustering to investigate further and to determine if we should remove any samples


# create a dataframe with samples as rows and genes as columns to input into hclust
countdata_forclustering <- t(countdata)
countdata_forclustering[1:5,1:10]

# run hierarchical clustering
dev.off()
hc <-hclust(dist(countdata_forclustering))
jpeg(filename= paste0(Output_Figures,"/",cur_date, "_Hierachical_clustering_outliers.jpeg"),width = 4500, height = 500)
plot(hc)
dev.off()

# from PCA and hierarchical clustering looks like 105_RNA-m105_Sal_24H_PLASMA_CATTTT_S39_L004 is an outlier



#################################################################################################
#################################################################################################
#### Removing sample(s) with potential QA/QC issues 
#################################################################################################
#################################################################################################


# First, pull all the IDs from the count data
keep_samples <- colnames(countdata)
# Note that in this case, this represents all 30 IDs

# Then, remove the one sample outlier from the "keep_samples" vectors
keep_samples <- keep_samples[! keep_samples %in% c("105_RNA-m105_Sal_24H_PLASMA_CATTTT_S39_L004")]
# Note that this represents all remaining 169 IDs

# Filter the countdata dataframe for these IDs
countdata <- countdata[, colnames(countdata) %in% keep_samples]
countdata[1:5, 1:10]

# Filter the coldata dataframe for these IDs
coldata <- coldata[coldata$Label_from_Core %in% keep_samples, ]


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
design <- model.matrix(~0 + Group + Batch, data=coldata)

#set up comparisons
#whatever is first in each statement is the numerator in the fold change estimate (treatment - control is like saying treatment/control)
#putting a plus sign in the formula looks at the additive affect of strain on the exposure effect; minus sign see the interaction of strain and exposure:
my_contrasts <- makeContrasts(RedOakSmolder_24h_Plasma = GroupRedOakSmolder_24h_Plasma - GroupSaline_24h_Plasma, 
                              PeatSmolder_24h_Plasma = GroupPeatSmolder_24h_Plasma - GroupSaline_24h_Plasma, 
                              RedOakFlame_24h_Plasma = GroupRedOakFlame_24h_Plasma - GroupSaline_24h_Plasma,
                              PeatFlame_24_Plasma = GroupPeatFlame_24h_Plasma - GroupSaline_24h_Plasma,
                              levels=design)


# Ensuring that the appropriate variable types are recognized within the DESeq2 algorithm (e.g., factor vs numeric)
coldata$Treatment <- as.factor(coldata$Treatment)
coldata$Group <- as.factor(coldata$Group)
coldata$Batch <- as.factor(coldata$Batch)



#########################################################
# Set up DESeq2 experiment
#########################################################
dds = DESeqDataSetFromMatrix(countData = countdata,
                             colData = coldata,
                             design = ~0+Group+Batch) #design matches the design above for specifying contrasts


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
