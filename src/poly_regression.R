#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library("ShortRead")

args <- commandArgs(trailingOnly = TRUE)

# Input FASTQ.
fl <- file.path(args[1])

# Output FASTQ
output_file <- file(args[2], "w")

# Set the min-max quality values.
min_max <- function(x) {
  if (x > 73)
    return(73)
  else if (x < 33)
    return(33)
  else
    return(x)
}

## Iterating over the FASTQ file.
f <- FastqStreamer(fl, 1)
while (length(fq <- yield(f))) {
  
  # Get the quality values.
  quals = as.numeric(as.vector(quality(fq)[[1]]))
  x = seq(1,length(quals))
  
  # Fit the polynomial function.
  fit = lm(quals ~ poly(x, strtoi(args[3]), raw=TRUE))
  
  # Predict the quality values using the above equation.
  predicted_quals = rawToChar(as.raw(unlist(lapply(round(predict(fit)), min_max))))
  #lines(x, round(predict(fit)))

  # Write out the fastq file to disk.
  # TODO: Replace with writeFastQ.
  writeLines(c(paste("@",as.vector(id(fq)), sep = ""), as.vector(sread(fq)), "+", predicted_quals), con=output_file)
}

close(output_file)
close(f)