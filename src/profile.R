#source("http://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library("ShortRead")

args <- commandArgs(trailingOnly = TRUE)

output_file <- file(args[2], "a")

training_set_size <- strtoi(args[3])

num_profiles <- strtoi(args[4])

# Taken from:
# http://stats.stackexchange.com/questions/78322/is-there-a-function-in-r-that-takes-the-centers-of-clusters-that-were-found-and
assign_to_cluster <- function(x, centers) {
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  max.col(-t(tmp))  # find index of min distance
}

# Read an entire fastq file
records <- readFastq(args[1])

# Take a training set, where we will run k-means on.
training_set = sample(records, size = training_set_size)

original_reads <- as(quality(records),'matrix')

# Run k-means.
mat = as(quality(training_set), 'matrix')
fit <- kmeans(mat, centers = num_profiles, iter.max = 600)

# Assign the original sequences to their closest cluster centers.
cluster_assignments = assign_to_cluster(original_reads, fit$centers)

# Write out each FASTQ record with its new profile quality values.
for (num in seq(1, length(records))) {
  predicted_quals = rawToChar(as.raw(round(fit$centers[cluster_assignments[num],] + 33)))
  writeLines(c(paste("@",as.vector(id(records[num])), sep = ""), as.vector(sread(records[num])), "+", as.vector(predicted_quals)), con=output_file)
}

close(output_file)