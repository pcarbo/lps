# Check in detail the topic proportions of the iLN samples. The last
# one (iLN_d2_20) is an outlier, and indeed this shows in the
# estimated topics.
library(ggplot2)
library(data.table)
library(fastTopics)
source("../code/lps_data.R")
set.seed(1)
cat("Loading count data.\n")
dat <- read_lps_data("../data/raw_read_counts.csv.gz")
samples <- dat$samples
counts  <- dat$counts
topic_colors <- c("darkblue","dodgerblue","darkorange","forestgreen",
                  "limegreen","tomato","darkred","olivedrab","magenta",
		  "darkmagenta","sienna","royalblue","lightskyblue",
                  "gold","red","cyan")
load("../output/fit-lps-k=16.RData")
fit <- poisson2multinom(fit)
i <- which(samples$tissue == "iLN")
fit1 <- select_loadings(fit,loadings = i)
structure_plot(fit1,loadings_order = 1:28,colors = topic_colors)
