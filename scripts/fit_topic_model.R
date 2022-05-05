# Fit a Poisson NMF model to the LPS RNA-seq data with k = 16 topics.
# These were the steps taken to load R and
# allocate computing resources for this analysis:
#
#   sinteractive -p broadwl -c 8 --mem=16G --time=72:00:00
#   model load R/3.5.1
#
library(tools)
library(data.table)
library(fastTopics)
source("../code/lps_data.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
cat("Loading count data.\n")
dat <- read_lps_data("../data/raw_read_counts.csv.gz")
samples <- dat$samples
counts  <- dat$counts

# Remove genes with very low (or no) expression.
j <- which(colSums(counts) > 20)
counts <- counts[,j]

# Fit a Poisson NMF model with k = 16 topics.
fit <- fit_poisson_nmf(counts,k = 16,init.method = "random",method = "em",
                       numiter = 20,control = list(nc = 8))
fit <- fit_poisson_nmf(counts,fit0 = fit,method = "scd",numiter = 180,
                       control = list(numiter = 4,nc = 8,extrapolate = TRUE))

# Perform DE analysis using the k = 16 topic model.
t0 <- proc.time()
de <- de_analysis(fit,counts,shrink.method = "ash",pseudocount = 0.1,
                  control = list(ns = 1e5,nc = 8))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Perform DE analysis using the topic model after merging some topics. 
fit_merged <- poisson2multinom(fit)
fit_merged <- merge_topics(fit_merged,c("k1","k5"))  # PBMC topics
fit_merged <- merge_topics(fit_merged,c("k2","k13")) # BM topics
fit_merged <- merge_topics(fit_merged,c("k6","k14")) # LI topics
t0 <- proc.time()
de_merged <- de_analysis(fit_merged,counts,shrink.method = "ash",
                         pseudocount = 0.1,control = list(ns = 1e5,nc = 8))
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save results to file.
save(list = c("fit","fit_merged","de","de_merged"),
     file = "fit-lps-k=16.RData")
resaveRdaFiles("fit-lps-k=16.RData")
