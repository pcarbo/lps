# A short script to compile a CSV file containing the results of the
# grade-of-membership differential expression (GoM DE) analysis based
# on the topic model with k = 16 topics.
library(fastTopics)
source("../code/de.R")

# Load the results of the topic modeling analysis.
load("../output/fit-lps-k=16.RData")

# Compile the DE results for all topics into a single table.
dat <- rbind(compile_de_table(de_merged),
             compile_de_table(de,c("k1","k2","k5","k6","k13","k14")))

# Filter out genes with lfsr >= 0.01.
dat <- subset(dat,lfsr < 0.01)

# Reorder the genes by topic, then by LFC.
topics <- c(paste0("k",1:16),c("k1+k5","k2+k13","k6+k14"))
dat$topic <- factor(dat$topic,topics)
rows <- with(dat,order(topic,-lfc))
dat  <- dat[rows,]

# Write the data frame to a CSV file.
dat <-
  transform(dat,
            topic = topic,
            lfc = format(round(lfc,digits = 3),trim = TRUE,scientific = FALSE),
            z = format(round(z,digits = 3),trim = TRUE,scientific = FALSE),
            lfsr = format(lfsr,digits = 3,trim = TRUE,scientific = TRUE))
write.csv(dat,"de_lps_k16.csv",quote = FALSE,row.names = FALSE)
