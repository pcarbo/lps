# TO DO: Explain here what this script is for, and how to use it.
#
#  sinteractive -p broadwl -c 8 --mem=8G --time=24:00:00
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

# Perform DE analysis using topic model.
# TO DO.

# Perform DE analysis using topic model with merged topics.
# TO DO.

# Save results to file.
save(list = "fit",file = "fit-lps-k=16.RData")
resaveRdaFiles("fit-lps-k=16.RData")
