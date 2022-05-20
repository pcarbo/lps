# A short script used to perform the gene set enrichment analysis
# using the GoM DE analysis results based on a topic model with k = 16
# topics.
library(Matrix)
library(tools)
library(susieR)
library(pathways)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# Load the gene sets.
data(gene_sets_mouse)
X <- gene_sets_mouse$gene_sets

# Load the results of the topic modeling analysis.
load("../output/fit-lps-k=16.RData")
topics <- c(paste0("k",1:16),c("k1+k5","k2+k13","k6+k14"))
Y <- cbind(de_merged$postmean,
           de$postmean[,c("k1","k2","k5","k6","k13","k14")])
Y <- Y[,topics]
           
# Remove gene sets in several MSigDB collections that clearly aren't
# relevant.
i <- which(!is.element(gene_sets_mouse$gene_set_info$database,
                       c("MSigDB-ARCHIVED","MSigDB-C1","MSigDB-C3",
                         "MSigDB-C4","MSigDB-C6")))
X <- X[,i]

# Align the gene-set data with the gene-wise statistics.
genes <- intersect(rownames(Y),gene_sets_mouse$gene_info$Symbol)
Y <- Y[genes,]
i <- match(genes,gene_sets_mouse$gene_info$Symbol)
X <- X[i,]
rownames(X) <- rownames(Y)
gene_sets_mouse$gene_info <- gene_sets_mouse$gene_info[i,]

# Next, remove gene sets with fewer than 10 genes and with more than
# 400 genes. Gene sets with a large number of genes are less likely to
# be interesting, and slow down the enrichment analysis, so they are
# removed.
i <- which(colSums(X) >= 10 & colSums(X) <= 400)
X <- X[,i]

# Perform a gene set enrichment analysis using susieR.
topics <- colnames(Y)
gsea <- vector("list",ncol(Y))
names(gsea) <- topics
t0 <- proc.time()
for (i in topics) {
  cat("topic",i,"\n")
  out <- susie(X,Y[,i],L = 10,intercept = TRUE,standardize = FALSE,
               estimate_residual_variance = TRUE,refine = FALSE,
               compute_univariate_zscore = FALSE,verbose = TRUE,
               min_abs_corr = 0)
  gsea[[i]] <- out[c("KL","lbf","sigma2","V","elbo","sets","pip","alpha","mu")]
}
t1 <- proc.time()
timing <- t1 - t0
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

# Save the results to file.
save(list = c("X","Y","gsea"),file = "gsea-lps-k=16.RData")
resaveRdaFiles("gsea-lps-k=16.RData")
