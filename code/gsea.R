# TO DO: Explain here what this function does, and how to use it.
compile_gsea_table <- function (s, gene_set_info) {

  # Initialize the output.
  out <- NULL
    
  # Get the credible sets.
  cs <- s$sets$cs

  # Add labels to some susie outputs.
  n                 <- length(s$lbf)
  names(s$lbf)      <- paste0("L",1:n)
  rownames(s$alpha) <- paste0("L",1:n)
  rownames(s$mu)    <- paste0("L",1:n)
  
  # Repeat for each CS.
  for (i in names(cs)) {
      
    # Get the variables (gene sets) included in the CS.
    j <- cs[[i]]
    n <- length(j)
    x <- data.frame(CS = rep(i,n),
                    lbf = rep(s$lbf[i],n),
                    stringsAsFactors = FALSE)
    x <- cbind(x,
               data.frame(pip  = s$alpha[i,j],
                          coef = s$mu[i,j]),
               gene_set_info[j,])
    out <- rbind(out,x)
  }

  rownames(out) <- NULL
  return(out)
}
