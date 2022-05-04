# TO DO: Explain here what this script is for, and how to use it.
library(data.table)
library(fastTopics)
library(ggplot2)
library(cowplot)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the count data.
dat     <- read_lps_data("../data/raw_read_counts.csv.gz")
samples <- dat$samples
genes   <- dat$genes
counts  <- dat$counts

# Remove genes with very low (or no) expression.
j <- which(colSums(counts) > 20)
counts <- counts[,j]

# Fit a Poisson NMF model with k = 16 topics.
fit <- fit_poisson_nmf(counts,k = 16,init.method = "random",method = "em",
                       numiter = 20,control = list(nc = 4))
fit <- fit_poisson_nmf(counts,fit0 = fit,method = "scd",numiter = 180,
                       control = list(numiter = 4,nc = 4,extrapolate = TRUE))

# Plot the improvement in the solution over time.
p1 <- plot_progress(fit,x = "timing",y = "loglik",colors = "black",
                    add.point.every = 10,e = 1e-4) +
  guides(color = "none",fill = "none",shape = "none",
         linetype = "none",size = "none")
p2 <- plot_progress(fit,x = "timing",y = "res",colors = "black",
                    add.point.every = 10,e = 1e-4) +
  guides(color = "none",fill = "none",shape = "none",
         linetype = "none",size = "none")
plot_grid(p1,p2)

# Create Structure plot.
set.seed(1)
pca_embed_method <- function (fit, ...)
  drop(pca_from_topics(fit,dims = 1))
p3 <- structure_plot(fit,grouping = samples$tissue,gap = 3,
                     colors = c("darkblue","dodgerblue","darkorange",
                                "forestgreen","limegreen","tomato","darkred",
                                "olivedrab","magenta","darkmagenta",
                                "sienna","royalblue","lightskyblue",
                                "gold","lightgray","cornflowerblue"),
                     topics = c(15,3,4,5,6,7,8,9,10,11,12,13,2,14,1,16),
                     embed_method = pca_embed_method)

# Perform DE analysis using topic model.
# TO DO.

# Save results to file.
# TO DO.
