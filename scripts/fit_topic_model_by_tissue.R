# TO DO: Explain here what this script is for, and how to use it.
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

# For each tissue, Fit a Poisson NMF model with k = 2 topics and k = 3
# topics.
tissues <- levels(samples$tissue)
n       <- length(tissues)
fits_k2 <- vector("list",n)
fits_k3 <- vector("list",n)
names(fits_k2) <- tissues
names(fits_k3) <- tissues
for (tissue in tissues) {
  cat(tissue,"\n")
  i <- which(samples$tissue == tissue)
  X <- counts[i,]

  # Remove genes with very low (or no) expression.
  j <- which(colSums(X) > 20)
  X <- X[,j]

  # Fit a Poisson NMF model with k = 2 topics.
  fit <- fit_poisson_nmf(X,k = 2,init.method = "random",method = "em",
                         numiter = 20,control = list(nc = 2),verbose = "none")
  fit <- fit_poisson_nmf(X,fit0 = fit,method = "scd",numiter = 180,
                         control = list(numiter = 4,nc = 2,extrapolate = TRUE),
                         verbose = "none")
  fits_k2[[tissue]] <- fit

  # Fit a Poisson NMF model with k = 3 topics.
  fit <- fit_poisson_nmf(X,k = 3,init.method = "random",method = "em",
                         numiter = 20,control = list(nc = 2),verbose = "none")
  fit <- fit_poisson_nmf(X,fit0 = fit,method = "scd",numiter = 180,
                         control = list(numiter = 4,nc = 2,extrapolate = TRUE),
                         verbose = "none")
  fits_k3[[tissue]] <- fit
}
