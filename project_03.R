#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ALDEx2")
#BiocManager::install("ANCOMBC")
#BiocManager::install("phyloseq")

# Load necessary libraries
library(phyloseq)
library(zCompositions)
# Load the compositions package
library(compositions)

# Imputation using mbImpute
# Handles zeros in the data by predicting missing values using a Bayesian framework
#install.packages("devtools")
#install.packages("glmnet")
#install.packages("devtools")
#install.packages("Matrix")
#library(devtools)
#install_github("lichen-lab/GMPR")
#install_github("ruochenj/mbImpute/mbImpute R package", force = TRUE)
# Load necessary libraries
library(GMPR)
library(mbImpute)

# Read the raw OTU table saved at genus level
feature_table_gen <- as.matrix(read.table("Genus_otu_table.txt", sep="\t", header = TRUE, row.names = 1))

# Read the metadata
metadata <- read.table("metadata_table.txt", sep="\t", header = TRUE, row.names = 1)

# CLR transformation with pseudocounts
# Adding a pseudocount to avoid issues with log of zero
feature_table_gen_pseudo <- feature_table_gen + 1
# Performing the clr transformation
clrtransform_pseudo <- apply(feature_table_gen_pseudo, 2, function(x) {
  log(x / exp(mean(log(x))))
})

# Optional: GMPR normalization
# GMPR is robust to zeros and does not assume a constant sum across samples
# Transpose the matrix so that samples are on rows, as expected by GMPR
transposed_feature_table <- as.data.frame(t(clrtransform_pseudo))
# Apply GMPR normalization
GMPR_factors <- GMPR(OTUmatrix = transposed_feature_table)
# Normalize the data
feature_table_gen_gmpr <- t(transposed_feature_table / GMPR_factors)

# mbImpute for imputation
# Handles zeros in the data by predicting missing values using a Bayesian framework
# mbImpute expects samples on rows... we need to transpose the count matrix
imp_count_mats <- mbImpute(condition = metadata$DiseaseState, otu_tab = t(feature_table_gen_gmpr), unnormalized = TRUE)
# mbImpute returns a list of three imputed OTU matrices. 
# imp_count_mat_lognorm: imputed normalized and log transformed matrix. 
# imp_count_mat_norm: imputed normalized count matrix with
# library size of each sample equal to 10^6. 
# imp_count_mat_origlibsize: imputed count matrix at the original library size.

# Selecting the imputed count matrix at the original library size
imp_count_mat <- t(imp_count_mats$imp_count_mat_origlibsize)

# Calculate the total number of zeros in the original dataset
tot_zeros <- sum(feature_table_gen == 0)

# Calculate the total number of zeros after imputation
tot_zeros_imp <- sum(imp_count_mat == 0) 

# Calculate and print the sparsity of the original dataset
sparsity <- round(tot_zeros / length(feature_table_gen) * 100, digits = 1)
paste("Sparsity =", sparsity, "%")

# Calculate and print the percentage of data imputed
perc_imputed <- round((tot_zeros - tot_zeros_imp) / length(feature_table_gen) * 100, digits = 1)
paste("Percentage imputed =", perc_imputed, "%")

# The choice to use clr transformation and mbImpute without normalization is based on the nature of the data and the goal to maintain the relative abundance information while addressing the issue of zeros in the dataset.

#imp_count_mats <- mbImpute(condition = label_samples$DiseaseState, otu_tab = t(feature_table_gen), unnormalized = T)

###################### 2 #############################

# Load necessary libraries

library(ALDEx2)
library(ANCOMBC)

# Assuming genus_otu_imputed and metadata are already loaded from Exercise 1

# Differential abundance analysis using Aldex2
# Aldex2 uses a Monte Carlo sampling technique to determine differences in abundance
aldex_results <- aldex(genus_otu_imputed, metadata$Condition, mc.samples = 128, denom = "all")

# Extracting differentially abundant taxa
aldex_signif <- aldex_results$we.eBH < 0.05

# Differential abundance analysis using Ancom
# Ancom assumes a log-normal distribution of the data
ancom_results <- ancom(phyloseq_object, "Condition")

# Extracting differentially abundant taxa
ancom_signif <- ancom_results$significance

# Compare the results
# Identify common and unique taxa
common_taxa <- intersect(names(aldex_signif), names(ancom_signif))
unique_aldex <- setdiff(names(aldex_signif), names(ancom_signif))
unique_ancom <- setdiff(names(ancom_signif), names(aldex_signif))

# Commenting on the results
# Number of DA taxa reported by each method and the common ones
num_aldex_da <- sum(aldex_signif)
num_ancom_da <- sum(ancom_signif)
num_common_da <- length(common_taxa)

# Properties/trends/characteristics of common taxa
# This requires further biological interpretation based on the taxa identified

# Output the results
list(
  num_aldex_da = num_aldex_da,
  num_ancom_da = num_ancom_da,
  num_common_da = num_common_da,
  common_taxa = common_taxa,
  unique_aldex = unique_aldex,
  unique_ancom = unique_ancom
)


#################### 3 #########################

# Load necessary libraries
library(vegan)
library(Rtsne)
library(umap)

# Assuming genus_otu_imputed, metadata, and DA taxa are already loaded from previous exercises

# Subset the data to include only DA taxa
genus_otu_da <- genus_otu_imputed[, DA_taxa]

# NMDS using the entire dataset
nmds_full <- metaMDS(genus_otu_imputed)
# NMDS using only DA taxa
nmds_da <- metaMDS(genus_otu_da)

# t-SNE using the entire dataset
tsne_full <- Rtsne(as.matrix(genus_otu_imputed))
# t-SNE using only DA taxa
tsne_da <- Rtsne(as.matrix(genus_otu_da))

# UMAP using the entire dataset
umap_full <- umap(as.matrix(genus_otu_imputed))
# UMAP using only DA taxa
umap_da <- umap(as.matrix(genus_otu_da))

# Plotting the results
par(mfrow=c(2,3))
plot(nmds_full, main="NMDS Full Dataset")
plot(nmds_da, main="NMDS DA Taxa")
plot(tsne_full$Y, main="t-SNE Full Dataset")
plot(tsne_da$Y, main="t-SNE DA Taxa")
plot(umap_full$layout, main="UMAP Full Dataset")
plot(umap_da$layout, main="UMAP DA Taxa")

# Commenting on the results
# Compare the plots visually and statistically (if applicable)
# Look for clustering patterns, overlap between groups, and separation of healthy donors and CDI patients

# Output the results
# Include any observations about trends or differences between the full dataset and the DA taxa
