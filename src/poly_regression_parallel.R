#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library("ShortRead")
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

# Input FASTQ.
fl <- file.path(args[1])

# Output FASTQ.
output_file <- file(args[2], "w")

# Polynomial degree.
degree <- args[3]

# Set the min-max quality values.
min_max <- function(x) {
  if (x > 73)
    return(73)
  else if (x < 33)
    return(33)
  else
    return(x)
}

poly_regression <- function(quals) {
  # Get the quality values. Use a hardcoded quality offset.
  quals = quals + 33
  x = seq(1,length(quals))
  
  # Fit the polynomial function.
  fit = lm(quals ~ poly(x, strtoi(degree)))
  
  # Predict the quality values using the above equation.
  predicted_quals = rawToChar(as.raw(unlist(lapply(round(predict(fit)), min_max))))
  return(predicted_quals)
}

# Read an entire fastq file
records <- readFastq(args[1])

original_reads <- as(quality(records),'matrix')

threads <- args[4]

cl <- makeCluster(strtoi(threads))
clusterExport(cl=cl, varlist=c("min_max", "degree"))

#ptm <- proc.time()
results <- parRapply(cl = cl, original_reads, poly_regression)

#print(proc.time() - ptm)
# poly_regression
# fptm <- proc.time()
# results2 <- apply(original_reads, MARGIN = 1, poly_regression)
# proc.time() - ptm

stopCluster(cl)

con  <- file(args[1], open = "r")
counter <- 1
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
  # Write the header...
  writeLines(oneLine, con=output_file)
  oneLine <- readLines(con, n = 1, warn = FALSE)
  
  # ... then the DNA sequence ...
  writeLines(oneLine, con=output_file)
  oneLine <- readLines(con, n = 1, warn = FALSE)
  
  # ... then the + ...
  writeLines(oneLine, con=output_file)
  oneLine <- readLines(con, n = 1, warn = FALSE)
  
  # ... and finally the new quality values.
  writeLines(results[[counter]], con=output_file)
  counter <- counter + 1
} 

close(con)