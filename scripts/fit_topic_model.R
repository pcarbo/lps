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
                     
