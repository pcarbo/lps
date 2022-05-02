library(data.table)
library(fastTopics)

# Load the count data.
counts <- fread("../data/raw_read_counts.csv.gz")
class(counts) <- "data.frame"
genes <- counts[,1]
counts <- t(as.matrix(counts[,-1]))
colnames(counts) <- genes
samples <- rownames(counts)
samples <- strsplit(samples,"_")
samples <- data.frame(tissue    = sapply(samples,"[[",1),
                      timepoint = sapply(samples,"[[",2),
                      mouse     = sapply(samples,"[[",3))

# Remove genes with very low (or no) expression.
j <- which(colSums(counts) > 20)
counts <- counts[,j]

# Fit a Poisson NMF model with k = 13 topics.
fit <- fit_poisson_nmf(counts,k = 13,init.method = "random",method = "em",
                       numiter = 20,control = list(nc = 2))
fit <- fit_poisson_nmf(counts,fit0 = fit,method = "scd",numiter = 180,
                       control = list(numiter = 4,nc = 2,extrapolate = TRUE))
