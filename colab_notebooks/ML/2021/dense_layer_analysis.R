#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(cowplot); library(ggplot2);

setwd("~/Desktop/ml_data/data_0527_9dil/0601_radial_log2_dense_layer_analysis/") # set to directory with CSV files containing imaging info
pref = '0601_radial_log2_mse_dense_'
plots <- list()
i = 1
for (n_layers in seq(from = 1, to = 10, by = 1)) {
  for (n_nodes in seq(from = 15, to = 50, by = 5)) {
    filename <- paste(pref, paste(toString(n_layers), paste('_', paste(toString(n_nodes), '.csv', sep = ""), sep = ""), sep = ""), sep = "")
    print(filename)
    k_df <- data.frame()
    k_df <- data.frame(read.csv(filename))
    plt <- ggplot(data = k_df) + geom_line(mapping = aes(X, X0), color = 'blue') + geom_line(mapping = aes(X, X0.1), color = 'orange') +
      labs(x = "epoch", y = "mse", title = paste("Layers: ", paste(toString(n_layers), paste(', Nodes: ', toString(n_nodes), sep = ""), sep = ""), sep = "")) +
      coord_cartesian(xlim = c(0, 300), ylim = c(0, 0.3))
    plots[[i]] <- plt
    i = i+1
  }
}
for (n_layers in seq(from = 1, to = 10, by = 1)) {
  print(plot_grid(plotlist = plots[((n_layers - 1) * 8 + 1):(n_layers * 8)]))
}

n_layers = 1
n_nodes = 15
filename <- paste('0525_mse_dense_', paste(toString(n_layers), paste('_', paste(toString(n_nodes), '.csv', sep = ""), sep = ""), sep = ""), sep = "")
print(filename)
k_df <- data.frame()
k_df <- data.frame(read.csv(filename))
ggplot(data = k_df) + geom_line(mapping = aes(X, X0), color = 'blue') + geom_line(mapping = aes(X, X0.1), color = 'orange') +
  labs(x = "epoch", y = "mse", title = paste("Layers: ", paste(toString(n_layers), paste(', Nodes: ', toString(n_nodes), sep = ""), sep = ""), sep = "")) +
  coord_cartesian(xlim = c(0, 300), ylim = c(0, 0.2))
plots[[i]] <- plt
i = i+1
