# ROC for analysis of SVM

#libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(pROC); library(MLmetrics); library(ggplot2);

setwd('~/Desktop/ml_data/data_0527')

# load in results of SVM
results <- data.frame(read.csv('results_0528_rf.csv'))
#response <- results$Actual
#predictor <- results$Prediction
response <- results$True
predictor <- results$Pred

roc <- pROC::roc(response, predictor)
plot(roc)
auc(roc)

MLmetrics::F1_Score(response, predictor)

proba <- data.frame(read.csv('proba_0528_rf.csv'))
ggplot(proba, aes(x=X0)) + geom_histogram(binwidth=0.05, color="black", fill = "white")

importances <- data.frame(read.csv('importances_0423.csv'))
feat_names <- data.frame(read.csv('~/Desktop/ml_data/feature_names.csv', header = FALSE))
feat_names <- feat_names[1:30,]
feat_imp <- cbind(feat_names, importances$X0)
colnames(feat_imp) <- c("feature", "importance")
feat_imp <- as.data.frame(feat_imp)
feat_imp <- feat_imp %>% dplyr::arrange(desc(importance))
plot(feat_imp$importance)

#perm_imp <- data.frame(read.csv('perm_importances_0415_100.csv'))
