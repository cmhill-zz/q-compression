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
  fit = lm(quals ~ poly(x, strtoi(degree), raw=TRUE))
  
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

# ptm <- proc.time()
results <- parApply(cl = cl, original_reads, MARGIN = 1, poly_regression)

# proc.time() - ptm
# poly_regression
# fptm <- proc.time()
# results2 <- apply(original_reads, MARGIN = 1, poly_regression)
# proc.time() - ptm

stopCluster(cl)

## Iterating over the FASTQ file again and print out the new quality values.
f <- FastqStreamer(fl, 1)
counter <- 1
while (length(fq <- yield(f))) {
  
  # Predict the quality values using the above equation.
  predicted_quals <- results[[counter]]
  
  # Write out the fastq file to disk.
  # TODO: Replace with writeFastQ.
  writeLines(c(paste("@",as.vector(id(fq)), sep = ""), as.vector(sread(fq)), "+", predicted_quals), con=output_file)
  
  counter <- counter + 1
}

close(output_file)
close(f)