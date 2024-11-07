# FIRST FUNCTION
MvAplot <- function(exprData, pdffilename, pcolor="black", lcolor="red") {
  # Open pdf file
  pdf(pdffilename)
  
  # Loop through each sample, excluding the first sample as reference
  for (i in 2:ncol(exprData)) {
    # Compute the M (log ratio) and A (average intensity) values
    # M = log-ratios=log2(r*1) - log2(r*2)
    # A = average =(log2(r*1) + log2(r*2)) / 2
    # CGL06 - Page 37
    M <- log2(exprData[, 1]) - log2(exprData[, i])
    A <- (log2(exprData[, 1]) + log2(exprData[, i])) / 2
    
    # Filter out NA and infinite values
    valid_idx <- is.finite(M) & is.finite(A)
    M <- M[valid_idx]
    A <- A[valid_idx]
    
    
    # Create the plot for each sample
    plot(A, M,
         main = paste("MvA Plot: Sample", colnames(exprData)[i], "vs Sample", colnames(exprData)[1]),
         xlab = "A (Average Intensity)",
         ylab = "M (Log Ratio)",
         col = pcolor, pch = 16, cex = 0.5)
    
    # Add a horizontal line at y = 0 with the specified line color
    abline(h = 0, col = lcolor, lwd = 2)
  }
  dev.off()
  print("PLOTTING IS DONE!")
}

# SECOND FUNCTION
TMMnorm <- function(exprData, annot, Mtrim=0.02, Atrim = c(0,8)) {
}

# THIRD FUNCTION
TMMnorm <- function(exprData, annot, Mtrim=0.02, Atrim = c(0,8)) {
  
}

# USAGE
rawData <- as.matrix(read.table("hw02/data/raw_count.txt", header = TRUE, row.names = 1, sep = "\t"))
# Add a small offset to avoid log2(0) issues
rawData_offset <- rawData + 1e-6

MvAplot(exprData = rawData, pdffilename = "hw02/MvAplots.pdf")
