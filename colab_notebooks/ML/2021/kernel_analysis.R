#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(cowplot);

setwd("~/Desktop/ml_data/data_0527_9dil/0601_radial_log2_kernel_analysis") # set to directory with CSV files containing imaging info

pref = '0601_radial_log2_mse_'

kernels = gsub(pattern = pref, replacement = "",x = gsub(pattern = ".csv", replacement = "", x = list.files(path = "./",pattern = "*.csv")))

plots <- list()
i = 1
for (k in kernels) {
  print(k)
  k_df <- data.frame(read.csv(paste(pref, paste(k, '.csv', sep = ""), sep = "")))
  k_reduced <- k_df[1:nrow(k_df),]
  plt <- ggplot(data = k_reduced) + geom_line(mapping = aes(X, X0), color = 'blue') + geom_line(mapping = aes(X, X0.1), color = 'orange') +
    labs(x = "epoch", y = "mse", title = k) +
    coord_cartesian(xlim = c(0, 200), ylim = c(0, 0.01))
  plots[[i]] <- plt
  i = i+1
}
plot_grid(plotlist = plots[1:4])

glorot_normal <- data.frame(read.csv('mse_glorot_normal.csv'))
glorot_normal_reduced <- glorot_normal[51:nrow(glorot_normal),]
ggplot(data = glorot_normal_reduced) + geom_line(mapping = aes(X, X0), color = 'blue') + geom_line(mapping = aes(X, X0.1), color = 'orange')
ggplot(data = k_reduced) + geom_line(mapping = aes(X, X0), color = 'blue') + geom_line(mapping = aes(X, X0.1), color = 'orange') +
  labs(x = "epoch", y = "mse", title = k)
