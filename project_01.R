analyze_SNP_data <- function(SNPfilepath, SNPannotationfilepath, indCTRL, MAFth=0.01, HWEalpha=0.01, alpha = 0.05) {
  # Read SNP annotation and SNP data files
  snp_annot <- read_snp_annotation(SNPannotationfilepath)
  snp_data <- read_snp_data(SNPfilepath)
  
  # Compute MAF for the SNP data
  maf_values <- compute_MAF(snp_data)
  
  # Compute HWE test p-values for the control group only
  hwe_p_values <- HWE_test(snp_data[, indCTRL])
  
  # Filter SNPs based on MAF and HWE p-values
  valid_snps <- names(maf_values)[maf_values >= MAFth & hwe_p_values >= HWEalpha]
  filtered_snp_data <- snp_data[valid_snps, , drop = FALSE]
  
  # Compute SNP-Phenotype association test p-values
  association_p_values <- SNP_association_test(filtered_snp_data, indCTRL)
  
  # Compute global p-values as the minimum of the three association p-values
  global_p_values <- apply(association_p_values, 1, min, na.rm = TRUE)
  
  # Perform multiple testing corrections
  num_tests <- length(global_p_values)
  
  # SIDAK correction
  sidak_threshold <- 1 - (1 - alpha)^(1 / num_tests)
  sidak_significant <- global_p_values <= sidak_threshold
  
  # Bonferroni correction
  bonferroni_threshold <- alpha / num_tests
  bonferroni_significant <- global_p_values <= bonferroni_threshold
  
  # Benjamini-Hochberg correction
  sorted_indices <- order(global_p_values)
  sorted_p_values <- global_p_values[sorted_indices]
  bh_q_values <- p.adjust(sorted_p_values, method = "BH")
  
  # Reorder bh_q_values to match original order of SNPs
  bh_q_values <- bh_q_values[order(sorted_indices)]
  
  # Associate SNPs with gene symbols, handling missing values
  gene_symbols <- snp_annot$Symbol[match(valid_snps, snp_annot$Pos)]
  gene_symbols[is.na(gene_symbols)] <- "NA"
  
  # Ensure consistency in vector lengths before creating the final data frame
  valid_snps <- intersect(valid_snps, rownames(association_p_values))
  gene_symbols <- gene_symbols[match(valid_snps, snp_annot$Pos)]
  gene_symbols[is.na(gene_symbols)] <- "NA"
  
  # Create the result dataframe
  result <- data.frame(
    Gene_Symbol = gene_symbols,
    MAF = maf_values[valid_snps],
    pvalHWE = hwe_p_values[valid_snps],
    pval_general = association_p_values[valid_snps, "pval_general"],
    pval_recessive = association_p_values[valid_snps, "pval_recessive"],
    pval_dominant = association_p_values[valid_snps, "pval_dominant"],
    pval = global_p_values[valid_snps],
    sidak = sidak_significant[valid_snps],
    bonferroni = bonferroni_significant[valid_snps],
    qval = bh_q_values[valid_snps],
    row.names = valid_snps
  )
  
  return(result)
}

# Example usage of the functions:

# Paths to the input files
data_file_path <- "./SNPdata.txt"  # Change with your path to data
annotation_file_path <- "./SNPAnnot.txt"  # Change with your path to annotation

# Example index of control subjects (modify according to your dataset)
indCTRL <- 1201:2000

# First test (MAFth=0.01, HWEalpha=0.01, alpha = 0.05)
result1 <- analyze_SNP_data(data_file_path, annotation_file_path, indCTRL, MAFth=0.01, HWEalpha=0.01, alpha = 0.05)
print(result1)

# Second test (MAFth=0.02, HWEalpha=0.02, alpha = 0.1)
result2 <- analyze_SNP_data(data_file_path, annotation_file_path, indCTRL, MAFth=0.02, HWEalpha=0.02, alpha = 0.1)
print(result2)
