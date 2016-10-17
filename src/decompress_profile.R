#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library("ShortRead")
library(parallel)

# Decompresses a binary-encoded profile created by profile_parallel.R.
# Rscript decompress_profile.R [1 compressed file] [2: num profiles] [3: read length] [4: num reads] [5: sequence fasta] [6: output name]

args <- commandArgs(trailingOnly = TRUE)
compressed_filename <- args[1]
num_profiles <- strtoi(args[2])
read_length <- strtoi(args[3])
num_reads <- strtoi(args[4])
sequence_filename <- args[5]
output_file <- file(args[6], 'w')

# If we need more than >256 profiles, we need more than one-byte of storage.
size <- 1
if (num_profiles > 256) {
  size <- 4
}

binary_input = file(compressed_filename, 'rb')
# Read in cluster center vectors.
vec <- readBin(binary_input, "integer", n = num_profiles * read_length, size = 1)
centers <- matrix(vec, nrow = num_profiles, ncol = read_length)
# Read in cluster assignments.
if (num_profiles > 256) {
  cluster_assignments <- readBin(binary_input, "integer", n = num_reads, size = size)
} else {
  cluster_assignments <- readBin(binary_input, "raw", n = num_reads, size = size)
}
cluster_assignments <- as.integer(cluster_assignments) + 1

con <- file(sequence_filename, open = 'r')
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
  writeLines(rawToChar(as.raw(round(centers[cluster_assignments[counter],]) + 33)), con=output_file)
  counter <- counter + 1
} 

close(output_file)
