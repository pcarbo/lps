# Read the count data from the CSV file (raw_read_counts.csv.gz).
read_lps_data <- function (file) {
  counts <- fread(file)
  class(counts) <- "data.frame"
  genes <- counts[,1]
  counts <- t(as.matrix(counts[,-1]))
  colnames(counts) <- genes
  samples <- rownames(counts)
  samples <- strsplit(samples,"_")
  samples <- data.frame(tissue    = sapply(samples,"[[",1),
                        timepoint = sapply(samples,"[[",2),
                        mouse     = sapply(samples,"[[",3))
  return(list(samples = samples,counts = counts))
}
