# fit <- poisson2multinom(fit)
pdat1 <- cbind(samples,data.frame(y = fit$L[,9]))
pdat2 <- cbind(samples,data.frame(y = fit$L[,1]))
pdat2 <- subset(pdat2,tissue == "PBMC")
pdat3 <- cbind(samples,data.frame(y = fit$L[,6]))
pdat3 <- subset(pdat3,tissue == "LI")
p1 <- ggplot(pdat1,aes(x = timepoint,y = y)) +
  geom_boxplot(width = 0.5) +
  labs(y = "topic 9",title = "all tissues")  +
  theme_cowplot(font_size = 12)
p2 <- ggplot(pdat2,aes(x = timepoint,y = y)) +
  geom_boxplot(width = 0.5) +
  labs(y = "topic 1",title = "PBMC") +
  theme_cowplot(font_size = 12)
p3 <- ggplot(pdat3,aes(x = timepoint,y = y)) +
  geom_boxplot(width = 0.5) +
  labs(y = "topic 6",title = "LI") +
  theme_cowplot(font_size = 12)
print(plot_grid(p1,p2,p3,nrow = 1,ncol = 3))
