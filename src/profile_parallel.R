#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library("ShortRead")
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
output_file <- file(args[2], "a")

training_set_size <- strtoi(args[3])

num_profiles <- strtoi(args[4])

# Taken from:
# http://stats.stackexchange.com/questions/78322/is-there-a-function-in-r-that-takes-the-centers-of-clusters-that-were-found-and
assign_to_cluster <- function(x) {
  # compute squared euclidean distance from each sample to each cluster center
   tmp <- apply(centers, 1, function(v) sum((x-v)^2))
   return(max.col(-t(tmp)))  # find index of min distance
}

# Read an entire fastq file
records <- readFastq(args[1])

# Take a training set, where we will run k-means on.
training_set = sample(records, size = training_set_size)

original_reads <- as(quality(records),'matrix')

# Run k-means.
mat = as(quality(training_set), 'matrix')
fit <- kmeans(mat, centers = num_profiles, iter.max = 600)
centers <- fit$centers

threads <- args[5]

cl <- makeCluster(strtoi(threads))
clusterExport(cl=cl, varlist=c("assign_to_cluster", "centers"))

#ptm <- proc.time()
# Assign the original sequences to their closest cluster centers.
cluster_assignments <- parRapply(cl = cl, original_reads, assign_to_cluster)

#print(proc.time() - ptm)
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
  #writeLines(results[[counter]], con=output_file)
  writeLines(rawToChar(as.raw(round(fit$centers[cluster_assignments[counter],] + 33))), con=output_file)
  counter <- counter + 1
} 

close(output_file)