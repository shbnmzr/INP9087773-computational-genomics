# Problem Description:
# Exercise 1: The occurrences observed on a sample of 1000 individuals shows aa genotype on 80
# subjects; AA on 500 subjects; Aa on 420 subjects.
# 1) Build a vector of numeric values with names corresponding to the genotype
# and values corresponding to occurrences
# 2) Calculates the minor allele frequency (MAF)


genotype_data <- c(500, 420, 80)
names(genotype_data) <- c("AA", "Aa", "aa")

compute_MAF <- function(occurrences) {
  if(!is.vector(occurrences)) {
    print("Not a valid input!")
    return()
  }
  
  total_count <- sum(occurrences)
  aa_count <- occurrences["aa"]
  Aa_count <- occurrences["Aa"]
  minor_allele_count <- (2 * aa_count) + Aa_count
  
  maf <- minor_allele_count / total_count
  return(maf)
}


# Test
compute_MAF(genotypes)