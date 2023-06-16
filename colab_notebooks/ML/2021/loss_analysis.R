#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(cowplot); library(ggplot2);

setwd("~/Desktop/ml_data/data_0527_9dil/0601_radial_log2_loss_analysis/") # set to directory with CSV files containing imaging info
pref = '0601_radial_log2_huber_mse_'
plots <- list()
i = 1
for (delta in seq(from = 0.1, to = 0.9, by = 0.1)) {
  filename <- paste(pref, paste(toString(delta), '.csv', sep = ""), sep = "")
  print(filename)
  k_df <- data.frame()
  k_df <- data.frame(read.csv(filename))
  plt <- ggplot(data = k_df) + geom_line(mapping = aes(X, X0), color = 'blue') + geom_line(mapping = aes(X, X0.1), color = 'orange') +
    labs(x = "epoch", y = "mse", title = paste("Delta: ", toString(delta), sep = "")) +
    coord_cartesian(xlim = c(0, 150), ylim = c(0, 0.2))
  plots[[i]] <- plt
  i = i+1
}
plot_grid(plotlist = plots[1:9])
