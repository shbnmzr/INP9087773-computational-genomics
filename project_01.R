compute_MAF <- function(SNPdata) {
  # Calculate the number of SNPs and subjects
  num_SNPs <- nrow(SNPdata)
  num_subjects <- ncol(SNPdata)
  
  # Print the number of SNPs and subjects
  cat("Number of SNPs:", num_SNPs, "\n")
  cat("Number of subjects:", num_subjects, "\n")
  
  # Initialize a vector to store MAF values
  MAF <- numeric(num_SNPs)
  
  # Calculate MAF for each SNP
  for (i in 1:num_SNPs) {
    # Count the number of each genotype
    genotypes <- factor(SNPdata[i, ], levels = c("0", "1", "2"))
    counts <- table(genotypes)
    
    # Ensure all genotypes are present in the counts
    counts <- counts + c("0" = 0, "1" = 0, "2" = 0)
    
    # Calculate frequencies of each allele
    # Assuming "0" is AA, "1" is Aa, and "2" is aa
    freq_AA <- counts["0"] / num_subjects
    freq_Aa <- counts["1"] / num_subjects
    freq_aa <- counts["2"] / num_subjects
    
    
    # Calculate MAF assuming 'aa' is the minor allele
    MAF[i] <- freq_aa + (freq_Aa / 2)
  }
  
  # Set the names of the MAF vector to match SNP IDs
  names(MAF) <- rownames(SNPdata)
  
  # Return the MAF vector
  return(MAF)
}



HWE_test <- function(SNPdata) {
  # Calculate the number of SNPs
  num_SNPs <- nrow(SNPdata)
  
  # Initialize a vector to store p-values
  p_values <- numeric(num_SNPs)
  
  # Loop through each SNP to perform HWE test
  for (i in 1:num_SNPs) {
    # Count the number of each genotype
    genotypes <- factor(SNPdata[i, ], levels = c("0", "1", "2"))
    counts <- table(genotypes)
    
    # Ensure all genotypes are present in the counts
    counts <- counts + c("0" = 0, "1" = 0, "2" = 0)
    
    
    # Observed genotype frequencies
    # Assuming "0" is AA, "1" is Aa, and "2" is aa
    obs_AA <- counts["0"]
    obs_Aa <- counts["1"]
    obs_aa <- counts["2"]
    
    # Total number of alleles
    total_alleles <- 2 * (obs_AA + obs_Aa + obs_aa)
    
    # Observed allele frequencies
    p_hat <- (2 * obs_AA + obs_Aa) / total_alleles
    q_hat <- (2 * obs_aa + obs_Aa) / total_alleles
    
    # Expected genotype frequencies under HWE
    exp_AA <- (p_hat^2) * (total_alleles / 2)
    exp_Aa <- 2 * p_hat * q_hat * (total_alleles / 2)
    exp_aa <- (q_hat^2) * (total_alleles / 2)
    
    # Chi-squared statistic
    chi_sq <- ((obs_AA - exp_AA)^2 / exp_AA) +
      ((obs_Aa - exp_Aa)^2 / exp_Aa) +
      ((obs_aa - exp_aa)^2 / exp_aa)
    
    # P-value calculation
    p_values[i] <- 1 - pchisq(chi_sq, df = 1)
  }
  
  # Set the names of the p-values vector to match SNP IDs
  names(p_values) <- rownames(SNPdata)
  
  # Find the SNP(s) with the lowest p-value
  min_p_value <- min(p_values)
  snps_with_min_p <- names(p_values)[p_values == min_p_value]
  
  # Print the SNP(s) with the lowest p-value
  cat("SNP(s) with the lowest p-value:", paste(snps_with_min_p, collapse = ", "), "\n")
  
  # Return the p-values vector
  return(p_values)
}



SNP_association_test <- function(SNPdata, indCTRL) {
  # Number of SNPs and subjects
  num_SNPs <- nrow(SNPdata)
  num_subjects <- ncol(SNPdata)
  
  # Initialize a matrix to store p-values
  p_values_matrix <- matrix(nrow = num_SNPs, ncol = 3)
  colnames(p_values_matrix) <- c("pval_general", "pval_recessive", "pval_dominant")
  
  # Calculate p-values for each SNP
  for (i in 1:num_SNPs) {
    # Split data into cases and controls
    controls <- SNPdata[i, indCTRL]
    cases <- SNPdata[i, -indCTRL]
    
    # General model
    genotypes <- factor(c(as.character(cases), as.character(controls)), levels = c("0", "1", "2"))
    contingency_table <- table(genotypes)
    observed <- c(contingency_table["0"], contingency_table["1"], contingency_table["2"])
    
    # Calculate allele frequencies
    p <- (2 * observed[1] + observed[2]) / (2 * sum(observed))
    q <- 1 - p
    
    # Expected genotype frequencies under HWE
    expected <- c(p^2, 2*p*q, q^2) * sum(observed)
    
    # Chi-squared statistic for general model
    chi_sq_general <- sum((observed - expected)^2 / expected)
    pval_general <- 1 - pchisq(chi_sq_general, df = 2)
    
    
    # Recessive model
    recessive_cases <- sum(cases == "2")
    recessive_controls <- sum(controls == "2")
    total_cases <- length(cases)
    total_controls <- length(controls)
    
    # Contingency table for recessive model
    contingency_table_recessive <- matrix(c(recessive_cases, total_cases - recessive_cases,
                                            recessive_controls, total_controls - recessive_controls),
                                          nrow = 2)
    
    # Chi-squared statistic for recessive model
    chi_sq_recessive <- chisq.test(contingency_table_recessive)$statistic
    pval_recessive <- chisq.test(contingency_table_recessive)$p.value
    
    
    # Dominant model
    dominant_cases <- sum(cases != "0")
    dominant_controls <- sum(controls != "0")
    
    # Contingency table for dominant model
    contingency_table_dominant <- matrix(c(dominant_cases, total_cases - dominant_cases,
                                           dominant_controls, total_controls - dominant_controls),
                                         nrow = 2)
    
    # Chi-squared statistic for dominant model
    chi_sq_dominant <- chisq.test(contingency_table_dominant)$statistic
    pval_dominant <- chisq.test(contingency_table_dominant)$p.value
    
    
    # Store p-values in the matrix
    p_values_matrix[i, ] <- c(pval_general, pval_recessive, pval_dominant)
  }
  
  # Set row names to SNP IDs
  rownames(p_values_matrix) <- rownames(SNPdata)
  
  # Plot boxplots of p-values
  par(mfrow = c(1, 3)) # Set layout for 3 boxplots in one panel
  boxplot(p_values_matrix[, "pval_general"], main = "General Model", ylab = "P-value")
  boxplot(p_values_matrix[, "pval_recessive"], main = "Recessive Model", ylab = "P-value")
  boxplot(p_values_matrix[, "pval_dominant"], main = "Dominant Model", ylab = "P-value")
  
  # Return the matrix of p-values
  return(p_values_matrix)
}

# Debugging steps inside the analyze_SNP_data() function



# And so on...


analyze_SNP_data <- function(SNPfilepath, SNPannotationfilepath, indCTRL, MAFth = 0.01, HWEalpha = 0.01, alpha = 0.05) {
  # Read the SNP data and annotation files
  SNPdata <- read.table(SNPfilepath, header = TRUE, row.names = 1)
  SNPannot <- read.table(SNPannotationfilepath, header = TRUE, sep = "\t", fill = TRUE)
  
  # Compute MAF, HWE p-values, and association test p-values
  MAF <- compute_MAF(SNPdata)
  pvalHWE <- HWE_test(SNPdata)
  
  # After computing MAF and HWE p-values
  print(head(MAF))
  print(head(pvalHWE))
  
  association_test_results <- SNP_association_test(SNPdata, indCTRL)
  
  # After computing association test p-values
  print(head(association_test_results))
  
  # Filter out SNPs based on MAF and HWE p-value
  filtered_SNPs <- which(MAF >= MAFth & pvalHWE >= HWEalpha)
  
  
  # After filtering SNPs
  print(filtered_SNPs)
  
  # Compute the global p-value
  pval <- apply(association_test_results[filtered_SNPs, 1:3], 1, min)
  
  # Perform multiple testing corrections
  sidak <- pval < (1 - (1 - alpha)^(1/length(pval)))
  bonferroni <- pval < (alpha / length(pval))
  qval <- p.adjust(pval, method = "BH")
  
  # Associate each SNP with its gene symbol
  gene_symbols <- SNPannot$Symbol[match(rownames(SNPdata)[filtered_SNPs], SNPannot$Pos)]
  
  # Create the output dataframe
  results_df <- data.frame(
    Gene_Symbol = gene_symbols,
    MAF = MAF[filtered_SNPs],
    pvalHWE = pvalHWE[filtered_SNPs],
    pval_general = association_test_results[filtered_SNPs, "pval_general"],
    pval_recessive = association_test_results[filtered_SNPs, "pval_recessive"],
    pval_dominant = association_test_results[filtered_SNPs, "pval_dominant"],
    pval = pval,
    sidak = sidak,
    bonferroni = bonferroni,
    qval = qval
  )
  
  # Set row names to SNP IDs
  rownames(results_df) <- rownames(SNPdata)[filtered_SNPs]
  
  return(results_df)
}



data_file_path = "SNPdata.txt" # change with your path to data
annotation_file_path = "SNPAnnot.txt" # change with your path to annotation

## First test (MAFth=0.01, HWEalpha=0.01, alpha = 0.05)

analyze_SNP_data(data_file_path, annotation_file_path, indCTRL=1201:2000, MAFth=0.01, HWEalpha=0.01, alpha = 0.05)
## NULL

## Second test (MAFth=0.02, HWEalpha=0.02, alpha = 0.1)

analyze_SNP_data(data_file_path, annotation_file_path, indCTRL=1201:2000, MAFth=0.02, HWEalpha=0.02, alpha = 0.1)
## NULL
