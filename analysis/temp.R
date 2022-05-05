pdat <- cbind(samples,
              data.frame(prop = poisson2multinom(fit)$L[,"k9"]))
ggplot(pdat,aes(x = timepoint,y = prop)) +
  geom_boxplot() +
  theme_cowplot()
