# TO DO: Explain here what this script does, and how to use it.
library(fastTopics)
source("../code/de.R")

# Load the results of the topic modeling analysis.
load("../output/fit-lps-k=16.RData")

# Compile the DE results for all topics into a single table.
dat <- compile_de_table(de_merged)

# Filter out genes with lfsr >= 0.01.
dat <- subset(dat,lfsr < 0.01)

# Write the data frame to a CSV file.
dat <-
  transform(dat,
            topic = topic,
            lfc = format(round(lfc,digits = 3),trim = TRUE,scientific = FALSE),
            z = format(round(z,digits = 3),trim = TRUE,scientific = FALSE),
            lfsr = format(lfsr,digits = 3,trim = TRUE,scientific = TRUE))
write.csv(dat,"de_table.csv",quote = FALSE,row.names = FALSE)
