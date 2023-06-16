# analyze results of regression

#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(ggplot2); library(cowplot);

setwd("~/Desktop/ml_data/HCR_0610/") # set to directory with CSV files containing imaging info

results_pred <- data.frame(read.csv('results_0610_507_D3_1_w1_limited.csv'))
results_true <- data.frame(read.csv('test_id_0610_507_D3_1_w1_limited.csv', header = FALSE))

results <- cbind(results_pred[,2:3], results_true)

colnames(results) <- c("Clu_Pred", "Aldob_Pred", "Clu_True", "Aldob_True")

results <- results %>% dplyr::mutate(cludif = Clu_True - Clu_Pred, entdif = Aldob_True - Aldob_Pred) %>%
                        dplyr::mutate(clusq = cludif * cludif, entsq = entdif * entdif) %>%
                        dplyr::mutate(magnitude = sqrt(clusq + entsq))
clumse <- (1 / nrow(results)) * sum(results$clusq)
entmse <- (1/ nrow(results)) * sum(results$entsq)
clumse
entmse

clu_r2 <- (cor(results$Clu_True, results$Clu_Pred))^2
aldob_r2 <- (cor(results$Aldob_True, results$Aldob_Pred))^2
clu_r2
aldob_r2

par(mfrow = c(2, 2))
plot(results$Clu_True, abs(results$cludif), xlab = 'Clu True', ylab = 'abs(residual)')
plot(results$Aldob_True, abs(results$entdif), xlab = "Aldob True", ylab = 'abs(residual)')
plot(results$Clu_True, results$Clu_Pred, xlab = "Clu True", ylab = "Clu Pred")
plot(results$Aldob_True, results$Aldob_Pred, xlab = "Aldob True", ylab = "Aldob Pred")

plots = list()
plots[[1]] <- ggplot(results, aes(Clu_True, Clu_Pred)) + geom_bin2d() + labs(x = 'Clu True', y = 'Clu Pred', title = "Clu") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plots[[2]] <- ggplot(results, aes(Aldob_True, Aldob_Pred)) + geom_bin2d() + labs(x = 'Aldob True', y = 'Aldob Pred', title = "Aldob") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plot_grid(plotlist = plots[1:2])

validation_results <- data.frame(read.csv('mse_0610_507_D3_1_w1_all.csv'))
ggplot(data = validation_results) + geom_line(mapping = aes(X, X0), color = 'blue') + geom_line(mapping = aes(X, X0.1), color = 'orange') +
  labs(x = "epoch", y = "mse") +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 0.05))




### Aldob only
results <- cbind(results_pred[,2], results_true)
colnames(results) <- c("Pred", "True")
results <- results %>% dplyr::mutate(entdif = Pred - True)
  dplyr::mutate(clusq = cludif * cludif, entsq = entdif * entdif) %>%
  dplyr::mutate(magnitude = sqrt(clusq + entsq))
clumse <- (1 / nrow(results)) * sum(results$clusq)
entmse <- (1/ nrow(results)) * sum(results$entsq)
clumse
entmse

par(mfrow = c(1, 1))
plot(results$Clu_True, abs(results$cludif), xlab = 'Clu True', ylab = 'abs(residual)')
plot(results$True, abs(results$entdif), xlab = "Aldob True", ylab = 'abs(residual)')
plot(results$Clu_True, results$Clu_Pred, xlab = "Clu True", ylab = "Clu Pred")
plot(results$True, results$Pred, xlab = "Aldob True", ylab = "Aldob Pred")

ggplot(results, aes(True, Pred)) + geom_bin2d() + labs(x = 'Aldob True', y = 'Aldob Pred', title = "Aldob") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

aldob_r2 <- (cor(results$True, results$Pred))^2


