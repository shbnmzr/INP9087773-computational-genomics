###################### 1 #############################

#Load the data
#read the raw otu table saved at genus level
feature_table_gen<-as.matrix(read.table("Genus_otu_table.txt",sep="\t"))

# Clr transformation with pseudocounts
feature_table_gen_pseudo<-feature_table_gen+1
np<-dim(feature_table_gen_pseudo)[2] #number of samples
nt<-dim(feature_table_gen_pseudo)[1] #number of taxa
clrtransform_pseudo<-feature_table_gen_pseudo #just to initialize clrtransform_pseudo
for (i in (1:np)) {
  den<-(prod(feature_table_gen_pseudo[,i]))^(1/nt) #geometric mean of column i
  clrtransform_pseudo[,i]<-log2(feature_table_gen_pseudo[,i]/den) #clr transformation of column i
}
head(clrtransform_pseudo)[1:5,1:10]

# Clr transformation without pseudocounts
np<-dim(feature_table_gen)[2] #number of samples
nt<-dim(feature_table_gen)[1] #number of taxa
clrtransform<-feature_table_gen #just to initialize clrtransform
for (i in (1:np)) {
  x<-feature_table_gen[,i]
  x[which(x==0)]<-NA
  den<-(prod(x,na.rm=TRUE)^(1/length(which(!is.na(x))))) #geometric mean of column i (excluding 0)
  clrtransform[,i]<-log2(x/den) #clr transformation of column i
  clrtransform[which(is.na(x)),i]<-0
}
head(feature_table_gen)[1:5,1:10]
head(clrtransform)[1:5,1:10]

#GMPR normalization
library(GMPR)
# GMPR expect samples on rows... we need to transpose the count matrix
GMPR_factors<- GMPR(OTUmatrix = as.data.frame(t(feature_table_gen)), min_ct = 2, intersect_no = 4) #see help for parameters meaning
feature_table_gen_gmpr<- t(t(feature_table_gen)/GMPR_factors)

GMPR_factors[1:10]
head(feature_table_gen)[1:5,1:10]
head(feature_table_gen_gmpr)[1:5,1:10]

#mbImpute (it takes around 8 minutes to run)
#it takes around 8 minutes to run
library(mbImpute)
label_samples<-read.table("metadata_table.txt",sep="\t",header=T)
# mbImpute expect samples on rows... we need to transpose the count matrix
imp_count_mats <- mbImpute(condition = label_samples$DiseaseState, otu_tab = t(feature_table_gen), unnormalized = T)
#mbImpute returns a list three imputed OTU matrices. 
#imp_count_mat_lognorm: imputed normalized and log transformed matrix. 
#imp_count_mat_norm : imputed normalized count matrix with
#library size of each sample equal to 10^6. 
#imp_count_mat_origlibsize: imputed countmatrix at the original library size.

imp_count_mat <-t(imp_count_mats[[3]])

#load("imp_count_mat_ls.RData")
#imp_count_mat <-t(imp_count_mat_ls[[3]])

#print original sparsity and % of imputed data
tot_zeros<- sum(feature_table_gen==0)
tot_zeros_imp<- sum(imp_count_mat==0) 

paste("Sparsity =", round(tot_zeros/length(feature_table_gen)*100, digits = 1))
paste("Perc imputed =",round((tot_zeros-tot_zeros_imp)/length(feature_table_gen)*100, digits = 1))

#PCA --- objects on rows so to project them on a lower dim space

PCA<-function(dati,condition){
  dati<-t(dati) # 
  N<-dim(dati)[1] #objects
  M<-dim(dati)[2] #genes (variables)
  S<-cov(dati)
  Eig<-eigen(S)
  lambda<-Eig[[1]] #eigenvalues
  PCs<-Eig[[2]] #eigenvectors (matrix V)
  varperc<-rep(0,M)
  for (i in (1:M)) varperc[i]<-sum(lambda[1:i])/sum(lambda)
  plot(varperc,type="b")
  Y<-dati%*%PCs # projection of the N objects in the new coordinates along the M PCs
  nmc<-names(table(condition))
  L<-length(nmc)
  plot(Y[,1],Y[,2]) # in this plot I am showing data in D=2 dim
  for (i in (2:L)) points(Y[which(condition==nmc[i]),1],Y[which(condition==nmc[i]),2],col=(i+1))
  #data can be reconstructed
  #rec_dati<-Y%*%t(PCs)
  return(list(Y,PCs,lambda))
}

#PCA of raw count data

resPCA<-PCA(feature_table_gen,condition=label_samples$DiseaseState)

#PCA of clr transformed data

resPCA<-PCA(clrtransform,condition=label_samples$DiseaseState)

#PCA of GMPR normalized data

resPCA<-PCA(feature_table_gen_gmpr,condition=label_samples$DiseaseState)

#PCA of mbImpute imputed data
resPCA<-PCA(imp_count_mat,condition=label_samples$DiseaseState)

#Impute the data and then clr transform to work on Euclidean space

# Clr transformation without pseudocounts on imputed data
np<-dim(imp_count_mat)[2] #number of samples
nt<-dim(imp_count_mat)[1] #number of taxa
impclr_count_mat<-imp_count_mat #just to initialize clrtransform
for (i in (1:np)) {
  x<-imp_count_mat[,i]
  x[which(x==0)]<-NA
  den<-(prod(x,na.rm=TRUE)^(1/length(which(!is.na(x))))) #geometric mean of column i (excluding 0)
  impclr_count_mat[,i]<-log2(x/den) #clr transformation of column i
  impclr_count_mat[which(is.na(x)),i]<-0
}
head(imp_count_mat)[1:5,1:10]
head(impclr_count_mat)[1:5,1:10]
resPCA<-PCA(impclr_count_mat,condition=label_samples$DiseaseState)

#Save workspace for the next hands-on

save.image("./preprocessing.RData")

###################### 2 #############################


# Load necessary libraries for DA analysis
library(ALDEx2)
library(ANCOMBC)
library(phyloseq)

# Assuming imp_count_mat and label_samples are already loaded from Exercise 1


# Differential abundance analysis using Aldex2
aldex.clr <- aldex.clr(imp_count_mat, label_samples$DiseaseState, mc.samples = 128, denom = "all")
aldex.effect <- aldex.effect(aldex.clr)
aldex.significant <- aldex.ttest(aldex.clr)
aldex.kw <- aldex.kw(aldex.clr)


# ♢
# rab.all - median clr value for all samples in the feature
# ♢
# rab.win.NS - median clr value for the NS group of samples
# ♢
# rab.win.S - median clr value for the S group of samples
# ♢
# rab.X1_BNS.q50 - median expression value of features in sample X1_BNS if include.item.summary=TRUE
# ♢
# dif.btw - median difference in clr values between S and NS groups
# ♢
# dif.win - median of the largest difference in clr values within S and NS groups
# ♢
# effect - median effect size: diff.btw / max(diff.win) for all instances
# ♢
# overlap - proportion of effect size distribution that overlaps 0 (i.e. no effect)
# ∗
# we.ep - Expected p-value of Welch’s t-test, a posterior predictive p-value
# ∗
# we.eBH - Expected Benjamini-Hochberg corrected p-value of Welch’s t test
# ∗
# wi.ep - Expected p-value of Wilcoxon rank test
# ∗
# wi.eBH - Expected Benjamini-Hochberg corrected p-value of Wilcoxon test

# Set the row names of imp_count_mat to match the IDs in label_samples
#flipped_imp_count_mat <- data.frame(t(imp_count_mat))
#rownames(imp_count_mat) <- label_samples$ID

# Create the OTU table for phyloseq
OTU <- otu_table(imp_count_mat, taxa_are_rows = TRUE)

# Create the sample data for phyloseq
sample_data <- sample_data(label_samples)

# Ensure that the sample names in the OTU table match the row names of the sample data
sample_names(OTU) <- rownames(sample_data)

# Now create the phyloseq object
physeq <- phyloseq(OTU, sample_data)

# Check if the phyloseq object is valid
physeq

# Run ANCOMBC2 for differential abundance analysis
ancombc2_res <- ancombc2(
  data = physeq,
  assay_name = "counts",
  tax_level = NULL, # Use NULL if you want to use the taxonomic level in the phyloseq object
  fix_formula = "DiseaseState", # Replace with the actual fixed effects formula
  rand_formula = NULL, # Use NULL if there are no random effects
  p_adj_method = "holm", # Choose an appropriate p-value adjustment method
  pseudo_sens = TRUE, # Set to TRUE to add a pseudocount for sensitivity
  prv_cut = 0.10, # Prevalence cut-off
  lib_cut = 1000, # Library size cut-off
  s0_perc = 0.05, # Threshold for the proportion of zero counts
  group = "DiseaseState", # Replace with the actual grouping variable
  struc_zero = FALSE, # Set to TRUE if structural zeros are expected
  neg_lb = FALSE, # Set to TRUE if negative binomial distribution is expected
  alpha = 0.05, # Significance level
  n_cl = 2, # Number of clusters for W-statistic calculation
  verbose = TRUE, # Set to TRUE for detailed output
  global = FALSE, # Set to TRUE for global test
  pairwise = TRUE, # Set to TRUE for pairwise comparisons
  dunnet = FALSE, # Set to TRUE for Dunnett's test
  trend = FALSE, # Set to TRUE for trend test
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE), # Iteration control for EM algorithm
  em_control = list(tol = 1e-5, max_iter = 100), # EM algorithm control
  lme_control = NULL, # Linear mixed-effects model control
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100), # Multiple testing control
  trend_control = NULL # Trend test control
)

# View the results
ancombc2_res

# merge into one output for convenience
aldex.all <- data.frame(aldex.significant,aldex.effect)

par(mfrow=c(1,3))
aldex.plot(aldex.all, type="MA", test="welch", main='MA plot')
aldex.plot(aldex.all, type="MW", test="welch", main='effect plot')
aldex.plot(aldex.all, type="volcano", test="welch", main='volcano plot')


# The left panel is the MA plot, the right is the MW (effect) plot. 
# In both plots red represents features called as differentially abundant with q <0.05; 
# grey are abundant, but not differentially abundant; black are rare, but not differentially abundant. 
# This function uses the combined output from the aldex.ttest and aldex.effect functions above



# Compare the results of the two methods
# Extracting significant features from both methods
aldex.sig.features <- rownames(aldex.significant)[aldex.significant$we.ep < 0.05]
ancombc2_res.sig.features <- rownames(ancombc2_res$res)[ancombc2_res$res$`p_(Intercept)` < 0.05]

# Remove "tax_" prefix
ancombc2_res.sig.features <- paste("tax_", ancombc2_res.sig.features, sep="")

# Identifying common features found by both methods
common.features <- intersect(aldex.sig.features, ancombc2_res.sig.features)

# Output the results
list(
  aldex_sig_features = aldex.sig.features,
  ancom_sig_features = ancombc2_res.sig.features,
  common_features = common.features
)

###################### 3 #############################
# Load necessary libraries for dimensionality reduction
library(vegan)
library(Rtsne)
library(uwot)

# Assuming imp_count_mat and label_samples are already loaded from Exercise 1

# Select DA taxa based on the results from Exercise 2
selected_taxa <- common.features

# Subset the data to include only DA taxa
da_taxa_data <- imp_count_mat[selected_taxa, ]

# NMDS using the DA taxa
nmds.da <- metaMDS(da_taxa_data)

# t-SNE using the DA taxa
tsne.da <- Rtsne(as.matrix(da_taxa_data))

# UMAP using the DA taxa
umap.da <- umap(as.matrix(da_taxa_data))

# Plotting the results
par(mfrow = c(1, 3))
plot(nmds.da, main = "NMDS on DA Taxa")
plot(tsne.da$Y, main = "t-SNE on DA Taxa")
plot(umap.da$layout, main = "UMAP on DA Taxa")

# Compare the plots with those obtained using all taxa
# You can run NMDS, t-SNE, and UMAP on the full dataset and compare visually
