# adapted from ml_celltype_dataprep.R

# used to look at feature space of nuclei and attribute these to particular celltype markers.
# prepares the data for use in machine learning.
# specifically prepares data for regression on features
# the difference is in how each cell is assigned a "percent" of a certain celltype
# allows for incorporation of data from multiple experiments--normalizes separately, combines to get training data

#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(tidyverse) ; library(ggplot2); library(cowplot);

setwd("~/Desktop/ml_data/HCR_0610/") # set to directory with CSV files containing imaging info

### Functions
getFrame = function(sample, local.radius=100, distal.radius = 500){
  sampledir = paste0(sample, "/")
  olfm4.df = data.frame(read.csv(paste0(sampledir,"1CY3_measurements.csv"))) # marker data for Cy3
  muc2.df = data.frame(read.csv(paste0(sampledir,"1YFP_measurements.csv"))) # marker data for Cy5
  chga.df = data.frame(read.csv(paste0(sampledir,"1CY5_measurements.csv"))) # marker data for GFP
  clu.df = data.frame(read.csv(paste0(sampledir,"2CY3_measurements.csv"))) # marker data for Cy3
  aldob.df = data.frame(read.csv(paste0(sampledir,"2YFP_measurements.csv"))) # marker data for Cy5
  lys1.df = data.frame(read.csv(paste0(sampledir,"2CY5_measurements.csv"))) # marker data for GFP

  nuclei = data.frame(read.csv(paste0(sampledir,"2_nuclei.feature.csv"))) # load in nulcear morphology data, may need to change number
  colnames(nuclei) <- c("subimage", "label", "area", "filled_area",
                        "bbox.0", "bbox.1", "bbox.2", "bbox.3", "centroid.0", "centroid.1",
                        "eccentricity", "solidity", "convex_area", "mean_intensity", "min_intensity", "max_intensity",
                        "orientation", "major_axis_length", "minor_axis_length", "perimeter", "extent", "intensity.image",
                        "near_contrast","near_dissim", "near_corr", "near_energy", "near_homog",
                        "mid_contrast","mid_dissim", "mid_corr", "mid_energy", "mid_homog",
                        "far_contrast","far_dissim", "far_corr", "far_energy", "far_homog",
                        "puncta_count", "puncta_area", "imvar")
  merge.df = merge(olfm4.df, merge(muc2.df, merge(chga.df, merge(clu.df, merge(aldob.df,lys1.df, by = c("subimage", "label")), by = c("subimage", "label")), by = c("subimage", "label")), by = c("subimage", "label")), by = c("subimage", "label")) # merge all channel data--add back other channels when needed
  #merge.df = merge(clu.df, aldob.df) # merge channel data--add back other channels when needed

  merge.final = merge(nuclei, merge.df, by = c("subimage", "label"))
  
  #merge.final <- data.frame(read.csv(paste0(sampledir, 'mergedMeasure.csv')))
  
  merge.final = merge.final[!(colnames(merge.final) %in% c("intensity.image")),]
  merge.final$sample = sample
  merge.final = merge.final[complete.cases(merge.final), ]

  merge.final$ratio = abs(merge.final$bbox.0 - merge.final$bbox.2)/abs(merge.final$bbox.1 - merge.final$bbox.3) # calculate aspect ratio
  # calculate minimum distance (centroid to centroid) between nuclei. also number of cells in a given radius.
  #sub.tab = data.frame(subset(merge.final, select = c("centroid.0","centroid.1")))
  # if (nrow(merge.final)<2){
  #   merge.final$mindist = 10**5
  #   merge.final$nearby.cells = 0
  #   merge.final$mid.cells = 0
  #   
  # }else{
  #   # look at distances between every pair of nuclei and take the minimum
  #   mindist = apply(data.matrix(as.matrix(dist(sub.tab))),1,FUN=function(x) {min(x[x > 0])})
  #   merge.final$mindist = mindist
  #   
  #   # treat every nucleus as a point (centroid), count all other nuclei within given radius
  #   nnframe = nn2(data = sub.tab,radius = local.radius, searchtype = 'radius', k = nrow(sub.tab))
  #   merge.final$nearby.cells = rowSums(data.matrix(nnframe$nn.idx[,-1]>0))
  #   
  #   nnframe = nn2(data = sub.tab,radius = local.radius, searchtype = 'radius', k = nrow(sub.tab))
  #   merge.final$mid.cells = rowSums(data.matrix(nnframe$nn.idx[,-1]>0))
  # }
  
  #dtab$mean_intensity = dtab$mean_intensity/median(dtab$mean_intensity)
  
  
  return(merge.final)
}

changeBounds=function(vec,log2_transform = T){
  if(log2_transform){
    vec = log2(vec)
  }
  maxV = max(vec)
  minV = min(vec)
  slope = 1/(maxV-minV)
  yInt = -1*slope*minV
  return(slope*vec+yInt)
}

### Execution

# get all image data names
sampleIDS = list.dirs(path = "./", full.names = FALSE, recursive = FALSE)
# construct large dataframe containing data for all images
all_data = c()
n_cells = c()

# use these later to get just the feature columns
featureCols = c("area", "filled_area","eccentricity", "solidity", "convex_area" #"mean_intensity", "min_intensity", "max_intensity",
                ,"orientation", "major_axis_length","minor_axis_length", "perimeter", "extent" 
                ,"near_contrast", "near_dissim","near_corr", "near_energy", "near_homog", "mid_contrast", "mid_dissim", "mid_corr", "mid_energy", "mid_homog", "far_contrast","far_dissim", "far_corr", "far_energy" 
                ,"puncta_count", "puncta_area",  "ratio")
idCols = c("subimage", "label")

# perform data transformation and normalization on individual samples/experiments, then combine into single dataframe
for (j in sampleIDS){
  print(j)
  final = getFrame(j)
  
  
  # filter out very large nuclei
  final = (final[complete.cases(final), ])
  final = final[final$area < 3000,]
  #final = final[final$mindist > 0,]
  
  #subdf<-subset.data.frame(final, select = c(idCols, featureCols))
  
  #feat_type = subdf
  feat_type <- final
  
  names(feat_type) <- names(feat_type) %>% gsub(pattern = "X1CY3", replacement = "Olfm4")  %>%
    gsub(pattern = "X1YFP", replacement = "Muc2")  %>%
    gsub(pattern = "X1CY5", replacement = "ChgA")  %>%
    gsub(pattern = "X2CY3", replacement = "Clu")  %>%
    gsub(pattern = "X2YFP", replacement = "Aldob")  %>%
    gsub(pattern = "X2CY5", replacement = "Lys1")  %>%
    gsub(pattern = "_mean_intensity", replacement = "")  %>%
    gsub(pattern = "X1YFP", replacement = "Muc2")  %>%
    gsub(pattern = "quantiles.0", replacement = "q10")  %>%
    gsub(pattern = "quantiles.1", replacement = "q25")  %>%
    gsub(pattern = "quantiles.2", replacement = "q50")  %>%
    gsub(pattern = "quantiles.3", replacement = "q75")  %>%
    gsub(pattern = "quantiles.4", replacement = "q90")
    
  # normalization on raw intensities
  feat_type$Olfm4_perc = ((feat_type$Olfm4 - min(feat_type$Olfm4))/(max(feat_type$Olfm4) - min(feat_type$Olfm4)))
  feat_type$Muc2_perc = ((feat_type$Muc2 - min(feat_type$Muc2))/(max(feat_type$Muc2) - min(feat_type$Muc2)))
  feat_type$ChgA_perc = ((feat_type$ChgA - min(feat_type$ChgA))/(max(feat_type$ChgA) - min(feat_type$ChgA)))
  feat_type$Clu_perc = ((feat_type$Clu - min(feat_type$Clu))/(max(feat_type$Clu) - min(feat_type$Clu)))
  feat_type$Aldob_perc = ((feat_type$Aldob - min(feat_type$Aldob))/(max(feat_type$Aldob) - min(feat_type$Aldob)))
  feat_type$Lys1_perc = ((feat_type$Lys1 - min(feat_type$Lys1))/(max(feat_type$Lys1) - min(feat_type$Lys1)))
  
  feat_type$distance = sqrt((feat_type$Clu_perc)^2 + (feat_type$Aldob_perc)^2)
  
  # log transform
  feat_type$Olfm4_log2trans = changeBounds(feat_type$Olfm4)
  feat_type$Muc2_log2trans = changeBounds(feat_type$Muc2)
  feat_type$ChgA_log2trans = changeBounds(feat_type$ChgA)
  feat_type$Clu_log2trans = changeBounds(feat_type$Clu)
  feat_type$Aldob_log2trans = changeBounds(feat_type$Aldob)
  feat_type$Lys1_log2trans = changeBounds(feat_type$Lys1)
  
  feat_type$log2distance = sqrt((feat_type$Clu_log2trans)^2 + (feat_type$Aldob_log2trans)^2) # change when using more channels
  
  feat_type$Clu_threshold = feat_type$Clu_q90/(feat_type$Clu_q50 + 1)
  feat_type$Aldob_threshold = feat_type$Aldob_q90/(feat_type$Aldob_q50 + 1)
  
  feat_type <- feat_type %>% dplyr::filter(Clu_threshold < 2 & Aldob_threshold < 2) # removes about 1000 cells
  feat_type <- feat_type %>% dplyr::filter(Clu < 1.1 * Clu_q50 & Clu > 0.9 * Clu_q50& Aldob < 1.1 * Aldob_q50 & Aldob > 0.9 * Aldob_q50) # removes about 5000 cells
  
  all_data <- rbind(all_data, feat_type)
  n_cells <- rbind(n_cells, as.data.frame(nrow(feat_type)))
 
}
plots = list()
plots[[1]] <- ggplot(all_data, aes(X1CY3_mean_intensity, X1CY3_quantiles.0)) + geom_bin2d() + geom_segment(aes(x = min(min(all_data$X1CY3_mean_intensity), min(all_data$X1CY3_quantiles.0)), y = min(min(all_data$X1CY3_mean_intensity), min(all_data$X1CY3_quantiles.0)), xend = max(max(all_data$X1CY3_mean_intensity), max(all_data$X1CY3_quantiles.0)), yend = max(max(all_data$X1CY3_mean_intensity), max(all_data$X1CY3_quantiles.0))), color = 'red')
plots[[2]] <- ggplot(all_data, aes(X1YFP_mean_intensity, X1YFP_quantiles.0)) + geom_bin2d() + geom_segment(aes(x = min(min(all_data$X1YFP_mean_intensity), min(all_data$X1YFP_quantiles.0)), y = min(min(all_data$X1YFP_mean_intensity), min(all_data$X1YFP_quantiles.0)), xend = max(max(all_data$X1YFP_mean_intensity), max(all_data$X1YFP_quantiles.0)), yend = max(max(all_data$X1YFP_mean_intensity), max(all_data$X1YFP_quantiles.0))), color = 'red')
plots[[3]] <- ggplot(all_data, aes(X1CY5_mean_intensity, X1CY5_quantiles.0)) + geom_bin2d() + geom_segment(aes(x = min(min(all_data$X1CY5_mean_intensity), min(all_data$X1CY5_quantiles.0)), y = min(min(all_data$X1CY5_mean_intensity), min(all_data$X1CY5_quantiles.0)), xend = max(max(all_data$X1CY5_mean_intensity), max(all_data$X1CY5_quantiles.0)), yend = max(max(all_data$X1CY5_mean_intensity), max(all_data$X1CY5_quantiles.0))), color = 'red')
plots[[4]] <- ggplot(all_data, aes(X2CY3_mean_intensity, X2CY3_quantiles.0)) + geom_bin2d() + geom_segment(aes(x = min(min(all_data$X2CY3_mean_intensity), min(all_data$X2CY3_quantiles.0)), y = min(min(all_data$X2CY3_mean_intensity), min(all_data$X2CY3_quantiles.0)), xend = max(max(all_data$X2CY3_mean_intensity), max(all_data$X2CY3_quantiles.0)), yend = max(max(all_data$X2CY3_mean_intensity), max(all_data$X2CY3_quantiles.0))), color = 'red')
plots[[5]] <- ggplot(all_data, aes(X2YFP_mean_intensity, X2YFP_quantiles.0)) + geom_bin2d() + geom_segment(aes(x = min(min(all_data$X2YFP_mean_intensity), min(all_data$X2YFP_quantiles.0)), y = min(min(all_data$X2YFP_mean_intensity), min(all_data$X2YFP_quantiles.0)), xend = max(max(all_data$X2YFP_mean_intensity), max(all_data$X2YFP_quantiles.0)), yend = max(max(all_data$X2YFP_mean_intensity), max(all_data$X2YFP_quantiles.0))), color = 'red')
plots[[6]] <- ggplot(all_data, aes(X2CY5_mean_intensity, X2CY5_quantiles.0)) + geom_bin2d() + geom_segment(aes(x = min(min(all_data$X2CY5_mean_intensity), min(all_data$X2CY5_quantiles.0)), y = min(min(all_data$X2CY5_mean_intensity), min(all_data$X2CY5_quantiles.0)), xend = max(max(all_data$X2CY5_mean_intensity), max(all_data$X2CY5_quantiles.0)), yend = max(max(all_data$X2CY5_mean_intensity), max(all_data$X2CY5_quantiles.0))), color = 'red')
plot_grid(plotlist=plots[4:5])

# heatmap of rawperc (only works for 2D)
ggplot(all_data, aes(Clu_perc, Aldob_perc)) + geom_bin2d()
ggplot(all_data, aes(Clu_log2trans, Aldob_log2trans)) + geom_bin2d()
ggplot(all_data, aes(Muc2_log2trans, Aldob_log2trans)) + geom_bin2d()
ggplot(all_data, aes(Lys1_log2trans, Aldob_log2trans)) + geom_bin2d()
ggplot(all_data, aes(ChgA_log2trans, Aldob_log2trans)) + geom_bin2d()

ggplot(all_data, aes(X2CY3_mean_intensity, X2CY3_quantiles.4)) + geom_bin2d() + coord_cartesian(xlim = c(0, 35000), ylim = c(0, 35000)) + geom_abline(aes(intercept = 0, slope = 1))
ggplot(all_data, aes(X2YFP_mean_intensity, X2YFP_quantiles.4)) + geom_bin2d() + coord_cartesian(xlim = c(0, 65000), ylim = c(0, 65000)) + geom_abline(aes(intercept = 0, slope = 1))

ggplot(feat_type %>% dplyr::filter(Clu < 1.1 * Clu_q50 & Clu > 0.9 * Clu_q50), aes(x = Clu, y = Clu - Clu_q50)) + geom_point()
ggplot(feat_type, aes(x = Clu, y = Clu - Clu_q50)) + geom_point()
ggplot(feat_type, aes(x = Aldob, y = Aldob - Aldob_q50)) + geom_point()
# radial quantiles (1/4 circles)
train_data = c()
test_data = c()
max = 7500
train = 6000
test = 1500
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  sampled_data <- all_data %>% dplyr::filter(distance > quantile(all_data$distance, probs = i) & distance <= quantile(all_data$distance, probs = i + 0.1)) %>%
    dplyr::sample_n(max, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):max,])
}

# radial with log2trans
# train_data = c()
# test_data = c()
# max = 3000
# train = 2700
# test = 300
# for (i in seq(from = 0, to = 0.9, by = 0.1)) {
#   sampled_data <- all_data %>% dplyr::filter(log2distance > quantile(all_data$log2distance, probs = i) & log2distance <= quantile(all_data$log2distance, probs = i + 0.1)) %>%
#     dplyr::sample_n(max, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data[1:train,])
#   test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }

changeBoundsQuantile=function(vec, log2_transform = F, upperPerc = 0.01, lowerPerc = 0.25) {
  if(log2_transform) {
    vec = log2(vec)
  }
  quant = quantile(vec, c(lowerPerc, 1 - upperPerc))
  qLower = quant[1]
  qUpper = quant[2]
  vec[vec < qLower] = qLower
  vec[vec > qUpper] = qUpper
  
  maxV = max(vec)
  minV = min(vec)
  slope = 1 / (maxV - minV)
  yInt = -1 * slope * minV
  return(slope * vec + yInt)
}
all_data$Clu_qtrans = changeBoundsQuantile(all_data$Clu)
all_data$Aldob_qtrans = changeBoundsQuantile(all_data$Aldob)
all_data <- all_data %>% dplyr::mutate(distance_qtrans = sqrt(Clu_qtrans^2 + Aldob_qtrans^2))
# radial, proportional to quantile size
sampleSubset = function(quantileNum, dfQuantile, dframe, split) {
  sampleRate = dfQuantile[dfQuantile$X1 == quantileNum,]$X2
  tmpFrame = dframe %>% dplyr::filter(quantile == quantileNum)
  total = floor(sampleRate * dim(tmpFrame)[1])
  train = as.integer(total * split)
  sampled_data = tmpFrame %>% dplyr::sample_n(total, replace = F)
  train_data  = sampled_data[1:train,]
  test_data = sampled_data[(train + 1):total,]
  return(list("train_data" = train_data, "test_data" = test_data))
}

subsample = function(vec = all_data, numQuantiles, maxSampleRate = 1, newTransform){
  if(newTransform) {
    vec$quantile <- ntile(vec$distance_qtrans, numQuantiles)
    quantList = quantile(vec$distance_qtrans, probs = seq(0, 1, 1/(numQuantiles-1)))
  } else {
    vec$quantile <- ntile(vec$distance, numQuantiles)
    quantList = quantile(vec$distance, probs = seq(0, 1, 1/(numQuantiles-1)))
  }
  trainingSampling = diff(c(0, quantList)) / max(c(0, diff(quantList))) * maxSampleRate
  dfQuant = data.frame(cbind(1:numQuantiles,as.numeric(trainingSampling)))
  train_data = c()
  test_data = c()
  for (i in 1:numQuantiles){
    data = sampleSubset(quantileNum = i, dfQuantile = dfQuant, dframe = vec, split = 0.8)
    train_data <- rbind(train_data, data$train_data)
    test_data <- rbind(test_data, data$test_data)
  }
  return(list("train_data" = train_data, "test_data" = test_data))
}

data = subsample(all_data, numQuantiles = 100, newTransform = T)
train_data <- data$train_data
test_data <- data$test_data
train_data <- subset(train_data, grepl('^\\d+$', train_data[,1]))
test_data <- subset(test_data, grepl('^\\d+$', test_data[,1]))
hist(train_data$distance, breaks = 100)

train_data = c()
test_data = c()
quants <- as.matrix(quantile(all_data$distance, probs = seq(0, 1, 0.1)))
quant_width <- c()
for (i in seq(0, 9, 1)) {
  quant_width[i] = quants[i + 1] - quants[i]
}
max_width <- max(quant_width)
sample_width <- quant_width / max_width
for (i in seq(0, 9, 1)) {
  range_data <- all_data %>% dplyr::filter(distance > quants[1, i / 10] & distance <= quants[1, i / 10 + 0.1])
  test = sample_width / 5
  sampled_data <- range_data %>% dplyr::sample_n(total, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[(test + 1):total,])
  test_data <- rbind(test_data, sampled_data[1:test,])
}
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  low_q = quantile(all_data$distance, probs = i)[[1]]
  upp_q = quantile(all_data$distance, probs = i + 0.1)[[1]]
  range_data <- all_data %>% dplyr::filter(distance > low_q & distance <= upp_q)
  scale_factor = 1 #reduce if too much data
  total = as.integer(nrow(range_data) * (upp_q - low_q) * scale_factor)
  print(total)
  test = total / 5
  print(test)
  sampled_data <- range_data %>% dplyr::sample_n(total, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[(test + 1):total,])
  test_data <- rbind(test_data, sampled_data[1:test,])
}

# use all data
total = nrow(all_data)
train = as.integer(0.8 * total)
shuffled_data <- all_data %>% sample_n(total, replace = FALSE)
train_data <- shuffled_data[1:train,]
test_data <- shuffled_data[(train+1):total,]


# plot test and train data
plots = list()
plots[[1]] <- ggplot(train_data, aes(Clu_perc, Aldob_perc)) + geom_bin2d() + labs(title = "Train Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plots[[2]] <- ggplot(test_data, aes(Clu_perc, Aldob_perc)) + geom_bin2d() + labs(title = "Test Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plot_grid(plotlist = plots[1:2])

plots = list()
plots[[1]] <- ggplot(train_data, aes(Clu_qtrans, Aldob_qtrans)) + geom_bin2d() + labs(title = "Train Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plots[[2]] <- ggplot(test_data, aes(Clu_qtrans, Aldob_qtrans)) + geom_bin2d() + labs(title = "Test Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plot_grid(plotlist = plots[1:2])
ggplot(all_data, aes(Clu_perc, Aldob_perc)) + geom_bin2d() + labs(title = "All Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

par(mfrow = c(1,2))
hist(train_data$Clu_qtrans)
hist(train_data$Aldob_qtrans)



featureColsLimited <- c("area", "eccentricity", "solidity", "perimeter", "major_axis_length")
range_data <- all_data %>% dplyr::filter((Clu_qtrans < 0.1 | Clu_qtrans > 0.9 | Aldob_qtrans < 0.1 | Aldob_qtrans > 0.9) & !(Clu_qtrans < 0.1 & Aldob_qtrans < 0.1))
lower_left <- all_data %>% dplyr::filter(Clu_qtrans < 0.1 & Aldob_qtrans < 0.1)
range_data <- rbind(range_data, lower_left %>% sample_n(as.integer(nrow(lower_left) * 0.1), replace = FALSE))
total = nrow(range_data)
train = as.integer(0.8 * total)
shuffled_data <- range_data %>% sample_n(total, replace = FALSE)
train_data <- shuffled_data[1:train,]
test_data <- shuffled_data[(train+1):total,]

range_data <- all_data %>% dplyr::filter(Clu_perc < 0.25 | Aldob_perc < 0.25 & !(Clu_perc < 0.25 & Aldob_perc < 0.25))
lower_left <- all_data %>% dplyr::filter(Clu_perc < 0.25 & Aldob_perc < 0.25)
range_data <- rbind(range_data, lower_left %>% sample_n(as.integer(nrow(lower_left) * 0.1), replace = FALSE))
total = nrow(range_data)
train = as.integer(0.8 * total)
shuffled_data <- range_data %>% sample_n(total, replace = FALSE)
train_data <- shuffled_data[1:train,]
test_data <- shuffled_data[(train+1):total,]


# write out the data
write.table(train_data %>% dplyr::select(all_of(featureColsLimited)), 
            file = "train_0610_507_D3_1_w1_limited.csv", row.names = F, col.names = F, sep = ",")
write.table(train_data %>% dplyr::select(Clu_perc, Aldob_perc), 
            file = "train_id_0610_507_D3_1_w1_limited.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(all_of(featureColsLimited)), 
            file = "test_0610_507_D3_1_w1_limited.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(Clu_perc, Aldob_perc), 
            file = "test_id_0610_507_D3_1_w1_limited.csv", row.names = F, col.names = F, sep = ",")
