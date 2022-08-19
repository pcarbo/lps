# TO DO: Explain here what this script does, and how to use it.
library(tools)
library(data.table)
library(fastTopics)
source("../code/lps_data.R")

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the sample meta data.
samples <-

# Load the count data.
dat <- read_lps_data("../data/cytokine_combo_2ndrun/counts_cytocombo.csv.gz")
counts  <- dat$counts
