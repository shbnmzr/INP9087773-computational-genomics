# Install edgeR if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("edgeR")

# Load the edgeR package
library(edgeR)




# Load the required data
exprData <- as.matrix(read.table("./raw_count.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = ""))
annot <- read.table("./gene_annot.txt", header = TRUE, sep = "\t", quote = "")

# Check the structure of the input data
str(exprData)
str(annot)

MvAplot <- function(exprData, pdffilename, pcolor = "black", lcolor = "red") {
  # Open a PDF device to save the plots
  pdf(pdffilename)
  
  # Get the number of samples
  num_samples <- ncol(exprData)
  
  # Loop over each sample (column) to create a plot against the first sample (column 1)
  for (i in 2:num_samples) {
    # Calculate M and A values
    A <- 0.5 * (log2(exprData[, 1] + 1) + log2(exprData[, i] + 1))  # Adding 1 to avoid log(0)
    M <- log2(exprData[, i] + 1) - log2(exprData[, 1] + 1)
    
    # Create the MvA plot
    plot(A, M, col = pcolor, pch = 20, 
         xlab = "A = 0.5 * (log2(Expression Sample 1 + 1) + log2(Expression Sample i + 1))", 
         ylab = "M = log2(Expression Sample i + 1) - log2(Expression Sample 1 + 1)",
         main = paste("MvA Plot: Sample 1 vs Sample", i))
    
    # Add horizontal line at y = 0
    abline(h = 0, col = lcolor, lwd = 2)
  }
  
  # Close the PDF device
  dev.off()
}

TMMnorm <- function(exprData, annot, Mtrim = 0.02, Atrim = c(0, 8)) {
  # Step 1: Scale the data by their sequencing depth and multiply by 10^6
  scaling_factors <- colSums(exprData)
  scaled_data <- t(t(exprData) / scaling_factors) * 1e6
  
  # Step 2: Calculate the scaling factors SF (with respect to sample 1)
  SF <- numeric(ncol(exprData))
  SF[1] <- 1  # Reference sample scaling factor is 1
  
  for (i in 2:ncol(exprData)) {
    M <- log2(scaled_data[, i] + 1) - log2(scaled_data[, 1] + 1)
    A <- 0.5 * (log2(scaled_data[, 1] + 1) + log2(scaled_data[, i] + 1))
    
    valid_indices <- which(A >= Atrim[1] & A <= Atrim[2])
    trimmed_M <- sort(M[valid_indices])
    trim_count <- ceiling(length(trimmed_M) * Mtrim)
    trimmed_M <- trimmed_M[(trim_count + 1):(length(trimmed_M) - trim_count)]
    
    SF[i] <- 2^mean(trimmed_M)
  }
  
  # Step 3: Normalize the data by their scaling factors SF with respect to sample 1
  normalized_data <- t(t(scaled_data) / SF)
  
  # Step 4: Scale the genes by their length and multiply by 10^3
  gene_lengths <- annot$Length[match(rownames(normalized_data), annot$Symbol)]
  normalized_data <- sweep(normalized_data, 1, gene_lengths, FUN = "/") * 1e3
  
  # Return the normalized matrix and scaling factors
  list(normalized_matrix = normalized_data, scaling_factors = SF)
}

DEbyEdgeR <- function(rawdata, groups, alpha = 0.05) {
  # Load necessary library
  library(edgeR)
  
  # Step 1: Create DGEList object
  dge <- DGEList(counts = rawdata, group = groups)
  
  # Step 2: Normalize the data using TMM normalization
  dge <- calcNormFactors(dge)
  
  # Step 3: Create the design matrix
  design <- model.matrix(~0 + groups)
  colnames(design) <- levels(factor(groups))
  
  # Step 4: Estimate dispersion
  dge <- estimateDisp(dge, design)
  
  # Step 5: Fit the model
  fit <- glmFit(dge, design)
  
  # Step 6: Perform likelihood ratio test (LRT)
  lrt <- glmLRT(fit, contrast = c(-1, 1))  # Compare all patients vs controls
  
  # Step 7: Extract p-values and log fold changes
  pvals <- lrt$table$PValue
  logFC <- lrt$table$logFC
  
  # Step 8: Calculate q-values using Benjamini-Hochberg procedure
  ranked_pvals <- rank(pvals, ties.method = "min")
  qvals <- pvals * length(pvals) / ranked_pvals
  qvals[qvals > 1] <- 1  # Ensure q-values are between 0 and 1
  
  # Step 9: Calculate the number of selected genes at significance level alpha
  selected_genes <- sum(pvals < alpha)
  
  # Step 10: Estimate false positives (FP), false negatives (FN), and FDR
  G <- nrow(rawdata)
  G0 <- 0.8 * G  # Assume 80% of genes are not differentially expressed
  FP <- alpha * G0
  FN <- (G - G0) * (1 - alpha)
  FDR <- FP / selected_genes
  
  # Step 11: Prepare the result matrix
  result_matrix <- data.frame(pval = pvals, qval = qvals, LFC = logFC)
  rownames(result_matrix) <- rownames(rawdata)
  
  # Return the result
  list(
    stats = c(
      selected_genes = selected_genes,
      FP = FP,
      FN = FN,
      FDR = FDR
    ),
    result_matrix = result_matrix
  )
}

DEbyWilcoxon <- function(rawdata, annotation, groups, alpha = 0.05, P0 = 0.8) {
  # Load necessary library
  library(edgeR)
  
  # Step 1: Normalize the input raw data using the TMM normalization
  norm_result <- TMMnorm(rawdata, annotation)
  norm_data <- norm_result$normalized_matrix
  
  # Step 2: Initialize variables to store results
  pvals <- rep(NA, nrow(norm_data))
  names(pvals) <- rownames(norm_data)
  
  # Step 3: Calculate the p-values using Wilcoxon test
  for (i in 1:nrow(norm_data)) {
    patient_data <- norm_data[i, groups == "patient"]
    control_data <- norm_data[i, groups == "CTRL"]
    
    # Ensure there are enough samples in both groups
    if (length(patient_data) > 1 && length(control_data) > 1) {
      test_result <- wilcox.test(patient_data, control_data)
      pvals[i] <- test_result$p.value
    }
  }
  
  # Step 4: Calculate q-values using Benjamini-Hochberg procedure
  qvals <- p.adjust(pvals, method = "BH")
  
  # Step 5: Calculate log fold change
  avg_patient <- rowMeans(norm_data[, groups == "patient"], na.rm = TRUE)
  avg_control <- rowMeans(norm_data[, groups == "CTRL"], na.rm = TRUE)
  LFC <- log2(avg_patient / avg_control)
  
  # Step 6: Calculate the number of DE genes
  num_DE_genes <- sum(qvals < alpha, na.rm = TRUE)
  
  # Step 7: Calculate the expected number of false positives (FP), false negatives (FN), and FDR
  G <- nrow(norm_data)
  G0 <- P0 * G
  FP <- alpha * G0
  FN <- G - G0 - num_DE_genes + FP
  FDR <- FP / max(1, num_DE_genes)
  
  # Step 8: Return the results
  result_matrix <- data.frame(pval = pvals, qval = qvals, LFC = LFC, row.names = rownames(norm_data))
  return(list(
    summary = c(num_DE_genes, FP, FN, FDR),
    result_matrix = result_matrix
  ))
}

# Install edgeR if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Update old packages
update.packages(ask = "a")

# Install the edgeR package using BiocManager
BiocManager::install("edgeR")

# Load the edgeR package
library(edgeR)

# Load the required data
raw_gene_expr <- as.matrix(read.table("./raw_count.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE, quote = ""))
gene_annotation <- read.table("./gene_annot.txt", header = TRUE, sep = "\t", quote = "")

# Vector of labels (first 19 columns are patient samples, the last 7 columns are control samples)
sample_groups <- factor(c(rep("patient", 19), rep("CTRL", 7)))

# DE using edgeR
edgeR_results <- DEbyEdgeR(rawdata = raw_gene_expr, groups = sample_groups, alpha = 0.05)

# DE using Wilcoxon
wilcox_results <- DEbyWilcoxon(rawdata = raw_gene_expr, annotation = gene_annotation, groups = sample_groups, alpha = 0.05, P0 = 0.8)

# Extract DE gene lists based on q-values < 0.05
edgeR_DE_genes <- rownames(edgeR_results$result_matrix)[edgeR_results$result_matrix$qval < 0.05]
wilcox_DE_genes <- rownames(wilcox_results$result_matrix)[wilcox_results$result_matrix$qval < 0.05]

# Number of DE genes detected by edgeR
num_DE_edgeR <- length(edgeR_DE_genes)

# Number of DE genes detected by Wilcoxon
num_DE_wilcox <- length(wilcox_DE_genes)

# Number of DE genes detected by both methods
common_DE_genes <- intersect(edgeR_DE_genes, wilcox_DE_genes)
num_common_DE_genes <- length(common_DE_genes)

# Print the results
cat("Number of DE genes detected by edgeR:", num_DE_edgeR, "\n")
cat("Number of DE genes detected by Wilcoxon:", num_DE_wilcox, "\n")
cat("Number of DE genes detected by both methods:", num_common_DE_genes, "\n")

# Check the trend in LFC between both methods
edgeR_LFC <- edgeR_results$result_matrix[common_DE_genes, "LFC"]
wilcox_LFC <- wilcox_results$result_matrix[common_DE_genes, "LFC"]

# Plot the LFC comparison
plot(edgeR_LFC, wilcox_LFC, xlab = "LFC by edgeR", ylab = "LFC by Wilcoxon", main = "Comparison of LFC between edgeR and Wilcoxon")
abline(0, 1, col = "red")

# Summary of trends in LFC
cat("Summary of LFC comparison:\n")
cat("Correlation between LFC values of common DE genes:", cor(edgeR_LFC, wilcox_LFC), "\n")
