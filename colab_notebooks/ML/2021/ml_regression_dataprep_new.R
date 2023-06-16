# adapted from ml_celltype_dataprep.R

# used to look at feature space of nuclei and attribute these to particular celltype markers.
# prepares the data for use in machine learning.
# specifically prepares data for regression on features
# the difference is in how each cell is assigned a "percent" of a certain celltype
# allows for incorporation of data from multiple experiments--normalizes separately, combines to get training data

#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(tidyverse) ; library(ggplot2); library(cowplot);

setwd("~/Desktop/ml_data/data_0527_9dil/") # set to directory with CSV files containing imaging info

### Functions
getFrame = function(sample, local.radius=100, distal.radius = 500){
  #a594.df = data.frame(read.csv(paste0(sample,"_a594_measurements.csv"))) # marker data for a594
  cy3.df = data.frame(read.csv(paste0(sample,"_cy3_measurements.csv"))) # marker data for Cy3
  cy5.df = data.frame(read.csv(paste0(sample,"_cy5_measurements.csv"))) # marker data for Cy5
  yfp.df = data.frame(read.csv(paste0(sample,"_yfp_measurements.csv"))) # marker data for GFP
  
  nuclei = data.frame(read.csv(paste0(sample,"_nuclei.feature.csv"))) # load in nulcear morphology data
  
  id = c(1:nrow(cy5.df))
  cy5.df = cbind(id, cy5.df)
  cy3.df = cbind(id, cy3.df)
  yfp.df = cbind(id, yfp.df)
  nuclei = cbind(id, nuclei)
  
  merge.df = merge(cy5.df, merge(cy3.df, yfp.df)) # merge all channel data--add back other channels when needed
  
  merge.final = merge(nuclei, merge.df, by = "id")
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
  print(slope)
  print(yInt)
  return(slope*vec+yInt)
}

### Execution

# get all image data names
sampleIDS = gsub(pattern = "_measurements.csv", replacement = "",x = gsub(pattern = "_yfp", replacement = "", x = list.files(path = "./",pattern = "*_yfp_measurements.csv")))
# construct large dataframe containing data for all images
all_data = c()
n_cells = c()
# perform data transformation and normalization on individual samples/experiments, then combine into single dataframe
for (j in sampleIDS){
  print(j)
  final = getFrame(j)
  
  
  # filter out very large nuclei
  final = (final[complete.cases(final), ])
  final = final[final$area < 3000,]
  #final = final[final$mindist > 0,]
  
  # these are the columns we use for PCA/UMAP
  goodCols = c("area", "filled_area","eccentricity", "solidity", "convex_area" #"mean_intensity", "min_intensity", "max_intensity",
               ,"orientation", "major_axis_length","minor_axis_length", "perimeter", "extent" 
               ,"near_contrast", "near_dissim","near_corr", "near_energy", "near_homog", "mid_contrast", "mid_dissim", "mid_corr", "mid_energy", "mid_homog", "far_contrast","far_dissim", "far_corr", "far_energy" 
               ,"puncta_count", "puncta_area",  "ratio")
  
  subdf<-subset.data.frame(final, select = goodCols)

  meanInt <- final %>% dplyr::select('cy5_mean_intensity', 'yfp_mean_intensity', 'cy3_mean_intensity') #add back other intensities as needed

  
  feat_type = subdf
  
  
  # original normalization
  feat_type$CY3 = (final$cy3_mean_intensity)#/min(feat_type$CY3_background))
  feat_type$CY5 = (final$cy5_mean_intensity)#/min(feat_type$CY5_background))
  feat_type$YFP = (final$yfp_mean_intensity)#/min(feat_type$YFP_background))
  feat_type$CY3minnorm = (feat_type$CY3/min(feat_type$CY3))
   feat_type$CY5minnorm = (feat_type$CY5/min(feat_type$CY5))
   feat_type$YFPminnorm = (feat_type$YFP/min(feat_type$YFP))
   feat_type$CY3norm = feat_type$CY3minnorm / (feat_type$CY3minnorm + feat_type$CY5minnorm + feat_type$YFPminnorm)
   feat_type$CY5norm = feat_type$CY5minnorm / (feat_type$CY3minnorm + feat_type$CY5minnorm + feat_type$YFPminnorm)
   feat_type$YFPnorm = feat_type$YFPminnorm / (feat_type$CY3minnorm + feat_type$CY5minnorm + feat_type$YFPminnorm)
  
  #plot(feat_type$CY5, feat_type$YFP)
  
  # re-normalization by percent
  # #feat_type$CY3perc = ((feat_type$CY3norm - min(feat_type$CY3norm))/(max(feat_type$CY3norm) - min(feat_type$CY3norm)))
  # feat_type$CY5perc = ((feat_type$CY5norm - min(feat_type$CY5norm))/(max(feat_type$CY5norm) - min(feat_type$CY5norm)))
  # feat_type$YFPperc = ((feat_type$YFPnorm - min(feat_type$YFPnorm))/(max(feat_type$YFPnorm) - min(feat_type$YFPnorm)))
  # #feat_type$CY3normperc = (feat_type$CY3perc / (feat_type$CY3perc + feat_type$CY5perc + feat_type$YFPperc))
  # feat_type$CY5normperc = (feat_type$CY5perc / (feat_type$CY3perc + feat_type$CY5perc + feat_type$YFPperc))
  # feat_type$YFPnormperc = (feat_type$YFPperc / (feat_type$CY3perc + feat_type$CY5perc + feat_type$YFPperc))
  
  # normalization on raw intensities
  feat_type$CY3rawperc = ((feat_type$CY3 - min(feat_type$CY3))/(max(feat_type$CY3) - min(feat_type$CY3)))
  feat_type$CY5rawperc = ((feat_type$CY5 - min(feat_type$CY5))/(max(feat_type$CY5) - min(feat_type$CY5)))
  feat_type$YFPrawperc = ((feat_type$YFP - min(feat_type$YFP))/(max(feat_type$YFP) - min(feat_type$YFP)))
  #feat_type$CY3rawnormperc = (feat_type$CY3rawperc / (feat_type$CY3rawperc + feat_type$CY5rawperc + feat_type$YFPrawperc))
  #feat_type$CY5rawnormperc = (feat_type$CY5rawperc / (feat_type$CY5rawperc + feat_type$YFPrawperc))
  #feat_type$YFPrawnormperc = (feat_type$YFPrawperc / (feat_type$CY5rawperc + feat_type$YFPrawperc))
  #plot(feat_type$CY5rawnormperc, feat_type$YFPrawnormperc)
  feat_type$CY3log2perc = log2(feat_type$CY3)
  feat_type$CY5log2perc = log2(feat_type$CY5)
  feat_type$YFPlog2perc = log2(feat_type$YFP)
  feat_type$CY3log10perc = log10(feat_type$CY3)
  feat_type$CY5log10perc = log10(feat_type$CY5)
  feat_type$YFPlog10perc = log10(feat_type$YFP)
  
  feat_type$distance = sqrt((feat_type$CY3rawperc)^2 + (feat_type$YFPrawperc)^2)
  
  feat_type$CY3log2trans = changeBounds(feat_type$CY3)
  feat_type$CY5log2trans = changeBounds(feat_type$CY5)
  feat_type$YFPlog2trans = changeBounds(feat_type$YFP)
  feat_type$log2distance = sqrt((feat_type$CY3log2trans)^2 + (feat_type$YFPlog2trans)^2)
  
  all_data <- rbind(all_data, feat_type)
  n_cells <- rbind(n_cells, as.data.frame(nrow(feat_type)))
  
  
}

# heatmap of rawperc (only works for 2D)
ggplot(all_data, aes(CY3rawperc,YFPrawperc)) + geom_bin2d()
ggplot(all_data, aes(CY3log2trans,YFPlog2trans)) + geom_bin2d()


par(mfrow = c(2, 3))
hist(all_data$CY3rawperc, ylim = range(0, 500), main = 'CY3 Raw')
hist(all_data$CY5rawperc, ylim = range(0, 500), main = 'CY5 Raw')
hist(all_data$YFPrawperc, ylim = range(0, 500), main = 'YFP Raw')

hist(all_data$CY3log2perc, breaks = 30, main = 'CY3 Log2')
hist(all_data$CY5log2perc, breaks = 30, main = 'CY5 Log2')
hist(all_data$YFPlog2perc, breaks = 30, main = 'YFP Log2')

hist(all_data$CY3log10perc, breaks = 30)
hist(all_data$CY5log10perc, breaks = 30)
hist(all_data$YFPlog10perc, breaks = 30)

hist(all_data$CY3rawperc, main = 'CY3 Raw')
hist(all_data$CY5rawperc, main = 'CY5 Raw')
hist(all_data$YFPrawperc, main = 'YFP Raw')

# choose one image for test, the rest for train
train_pool <- all_data[25033:nrow(all_data),]
test_pool <- all_data[1:25032,]
# train_data <- train_pool %>% dplyr::sample_n(1000, replace = FALSE)
# test_data <- test_pool %>% dplyr::sample_n(100, replace = FALSE)

# prepare data for ml
# first need to evaluate data distribution
# sample evenly across intensities -- will need to adjust the number
### 1D sampling
# max = 32
# train = 24
# test = 8
# sampled_data <- all_data %>% dplyr::filter(CY5rawperc >= 0.8 & CY5rawperc <= 1) %>%
#                   dplyr::sample_n(max, replace = FALSE)
# train_data <- sampled_data[1:train,]
# test_data <- sampled_data[(train + 1):max,]
# for (i in seq(from = 0, to = 0.6, by = 0.2)) {
#   sampled_data <- all_data %>% dplyr::filter(CY5rawperc >= i & CY5rawperc < i + 0.2) %>%
#     dplyr::sample_n(max, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data[1:train,])
#   test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc >= 0.8 & YFPrawperc <= 1) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# for (i in seq(from = 0, to = 0.6, by = 0.2)) {
#   sampled_data <- all_data %>% dplyr::filter(YFPrawperc >= i & YFPrawperc < i + 0.2) %>%
#     dplyr::sample_n(max, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data[1:train,])
#   test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }


### random sample
# sampled_data <- all_data %>% dplyr::sample_n(2200, replace = FALSE)
# train_data <- sampled_data[1:2000,]
# test_data <- sampled_data[2001:2200,]

###2D sampling (even binning)
# max = 20
# train = 15
# test = 5
# sampled_data <- all_data %>% dplyr::filter(CY3rawperc >= 0.6 & CY3rawperc <= 1 & YFPrawperc >= 0 & YFPrawperc < 0.2) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- sampled_data[1:train,]
# test_data <- sampled_data[(train + 1):max,]
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc >= 0.6 & YFPrawperc <= 1 & CY3rawperc >= 0 & CY3rawperc < 0.2) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# for (i in seq(from = 0.0, to = 0.4, by = 0.2)) {
#   sampled_data <- all_data %>% dplyr::filter(CY3rawperc >= i & CY3rawperc <= i + 0.2 & YFPrawperc < 0.2) %>%
#     dplyr::sample_n(max, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data[1:train,])
#   test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }
# for (i in seq(from = 0.2, to = 0.4, by = 0.2)) {
#   sampled_data <- all_data %>% dplyr::filter(YFPrawperc >= i & YFPrawperc <= i + 0.2 & CY5rawperc < 0.2) %>%
#     dplyr::sample_n(max, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data[1:train,])
#   test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }
# max = 8
# train = 6
# test = 2
# for (i in seq(from = 0.2, to = 0.8, by = 0.2)) {
#   sampled_data <- all_data %>% dplyr::filter(CY5rawperc >= i & CY5rawperc < i + 0.2 & YFPrawperc >= 0.2 & YFPrawperc < 0.4) %>%
#     dplyr::sample_n(max, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data[1:train,])
#   test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }
# for (i in seq(from = 0.4, to = 0.8, by = 0.2)) {
#   sampled_data <- all_data %>% dplyr::filter(YFPrawperc >= i & YFPrawperc < i + 0.2 & CY5rawperc >= 0.2 & CY5rawperc < 0.4) %>%
#     dplyr::sample_n(max, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data[1:train,])
#   test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }
# 

# 
# ### with 2D concentrated sampling (binning)
# train_data = c()
# test_data = c()
# max = 20
# train = 15
# test = 5
# # size 1 (1/16)
# for (i in seq(from = 0.0, to = 0.15, by = 0.05)) {
#   for (j in seq(from = 0.0, to = 0.15, by = 0.05)) {
#     if (i >= 0.1 & j >= 0.1) {
#       break
#     }
#     sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.05 & CY3rawperc > j & CY3rawperc <= j + 0.05) %>%
#       dplyr::sample_n(max, replace = FALSE)
#     train_data <- rbind(train_data, sampled_data[1:train,])
#     test_data <- rbind(test_data, sampled_data[(train + 1):max,])
#   }
# }
# # size 2 (1/4)
# i = 0.1
# j = 0.1
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.1 & CY3rawperc > j & CY3rawperc <= j + 0.1) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# 
# i = 0.2
# for (j in c(0, 0.1)) {
#     print(i)
#     print(j)
#     sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.1 & CY3rawperc > j & CY3rawperc <= j + 0.1) %>%
#       dplyr::sample_n(max, replace = FALSE)
#     train_data <- rbind(train_data, sampled_data[1:train,])
#     test_data <- rbind(test_data, sampled_data[(train + 1):max,])
#     sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.1 & CY3rawperc > i & CY3rawperc <= i + 0.1) %>%
#       dplyr::sample_n(max, replace = FALSE)
#     train_data <- rbind(train_data, sampled_data[1:train,])
#     test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# }
# 
# # # size 3 (1)
# # i = 0.4
# # j = 0.0
# # sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.2 & CY3rawperc > j & CY3rawperc <= j + 0.2) %>%
# #   dplyr::sample_n(max, replace = FALSE)
# # train_data <- rbind(train_data, sampled_data[1:train,])
# # test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# # sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.2 & CY3rawperc > i & CY3rawperc <= i + 0.2) %>%
# #   dplyr::sample_n(max, replace = FALSE)
# # train_data <- rbind(train_data, sampled_data[1:train,])
# # test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# # i = 0.4
# # j = 0.2
# # sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.2 & CY3rawperc > j & CY3rawperc <= j + 0.2) %>%
# #   dplyr::sample_n(max, replace = FALSE)
# # train_data <- rbind(train_data, sampled_data[1:train,])
# # test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# # sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.2 & CY3rawperc > i & CY3rawperc <= i + 0.2) %>%
# #   dplyr::sample_n(max, replace = FALSE)
# # train_data <- rbind(train_data, sampled_data[1:train,])
# # test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# # i = 0.2
# # j = 0.2
# # sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.2 & CY3rawperc > j & CY3rawperc <= j + 0.2) %>%
# #   dplyr::sample_n(max, replace = FALSE)
# # train_data <- rbind(train_data, sampled_data[1:train,])
# # test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# 
# # size 3.5 (1)
# i = 0.3
# j = 0.0
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.3 & CY3rawperc > j & CY3rawperc <= j + 0.3) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.3 & CY3rawperc > i & CY3rawperc <= i + 0.3) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# 
# # size 4 (4)
# i = 0.6
# j = 0.0
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.4 & CY3rawperc > j & CY3rawperc <= j + 0.4) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.4 & CY3rawperc > i & CY3rawperc <= i + 0.4) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# # size 5 (9)
# # sampled_data <- all_data %>% dplyr::filter(YFPrawperc > 0.4 & YFPrawperc <= 1 & CY3rawperc > 0.4 & CY3rawperc <= 1) %>%
# #   dplyr::sample_n(max, replace = FALSE)
# # train_data <- rbind(train_data, sampled_data[1:train,])
# # test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# 
# # undersampling bottom
# train_data <- c()
# test_data <- c()
# max = 20
# train = 15
# test = 5
# i = 0
# j = 0
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.3 & CY3rawperc > j & CY3rawperc <= j + 0.3) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# j = 0.3
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.3 & CY3rawperc > j & CY3rawperc <= j + 0.3) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.3 & CY3rawperc > i & CY3rawperc <= i + 0.3) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# i = 0.6
# j = 0.0
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.4 & CY3rawperc > j & CY3rawperc <= j + 0.4) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.4 & CY3rawperc > i & CY3rawperc <= i + 0.4) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# 
# # oversampling middle
# train_data <- c()
# test_data <- c()
# max = 20
# train = 15
# test = 5
# # i = 0
# # j = 0
# # sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.3 & CY3rawperc > j & CY3rawperc <= j + 0.3) %>%
# #   dplyr::sample_n(max, replace = FALSE)
# # train_data <- rbind(train_data, sampled_data[1:train,])
# # test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# i = 0.6
# j = 0.0
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.4 & CY3rawperc > j & CY3rawperc <= j + 0.4) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.4 & CY3rawperc > i & CY3rawperc <= i + 0.4) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# i = 0
# j = 0.3
# max = 100
# train = 75
# test = 25
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > i & YFPrawperc <= i + 0.3 & CY3rawperc > j & CY3rawperc <= j + 0.3) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])
# sampled_data <- all_data %>% dplyr::filter(YFPrawperc > j & YFPrawperc <= j + 0.3 & CY3rawperc > i & CY3rawperc <= i + 0.3) %>%
#   dplyr::sample_n(max, replace = FALSE)
# train_data <- rbind(train_data, sampled_data[1:train,])
# test_data <- rbind(test_data, sampled_data[(train + 1):max,])

# radial quantiles (1/4 circles)
train_data = c()
test_data = c()
max = 3000
train = 2700
test = 300
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  sampled_data <- all_data %>% dplyr::filter(distance > quantile(all_data$distance, probs = i) & distance <= quantile(all_data$distance, probs = i + 0.1)) %>%
    dplyr::sample_n(max, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):max,])
}
# radial with log2trans
train_data = c()
test_data = c()
max = 3000
train = 2700
test = 300
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  sampled_data <- all_data %>% dplyr::filter(log2distance > quantile(all_data$log2distance, probs = i) & log2distance <= quantile(all_data$log2distance, probs = i + 0.1)) %>%
    dplyr::sample_n(max, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):max,])
}

# bisected radial quantiles (1/8 circles)
train_data = c()
test_data = c()
max = 1232
train = 1109
test = 123
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  upper <- all_data %>% dplyr::filter(YFPrawperc > CY3rawperc)
  sampled_data <- upper %>% dplyr::filter(distance > quantile(upper$distance, probs = i) & distance <= quantile(upper$distance, probs = i + 0.1)) %>%
    dplyr::sample_n(max, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):max,])
}
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  lower <- all_data %>% dplyr::filter(CY3rawperc > YFPrawperc)
  sampled_data <- lower %>% dplyr::filter(distance > quantile(lower$distance, probs = i) & distance <= quantile(lower$distance, probs = i + 0.1)) %>%
    dplyr::sample_n(max, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):max,])
}

## radial but separate test and train images
# radial quantiles (1/4 circles)
# train_data = c()
# test_data = c()
# train = 2000
# test = 200
# for (i in seq(from = 0, to = 0.9, by = 0.1)) {
#   sampled_data <- train_pool %>% dplyr::filter(distance > quantile(train_pool$distance, probs = i) & distance <= quantile(train_pool$distance, probs = i + 0.1)) %>%
#     dplyr::sample_n(train, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data)
#   sampled_data <- test_pool %>% dplyr::filter(distance > quantile(test_pool$distance, probs = i) & distance <= quantile(test_pool$distance, probs = i + 0.1)) %>%
#     dplyr::sample_n(test, replace = FALSE)
#   test_data <- rbind(test_data, sampled_data)
# }
# 
# # bisected radial quantiles (1/8 circles)
# train_data = c()
# test_data = c()
# train = 600
# test = 100
# for (i in seq(from = 0, to = 0.9, by = 0.1)) {
#   upper <- train_pool %>% dplyr::filter(YFPrawperc > CY3rawperc)
#   sampled_data <- upper %>% dplyr::filter(distance > quantile(upper$distance, probs = i) & distance <= quantile(upper$distance, probs = i + 0.1)) %>%
#     dplyr::sample_n(train, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data)
#   upper <- test_pool %>% dplyr::filter(YFPrawperc > CY3rawperc)
#   sampled_data <- upper %>% dplyr::filter(distance > quantile(upper$distance, probs = i) & distance <= quantile(upper$distance, probs = i + 0.1)) %>%
#     dplyr::sample_n(test, replace = FALSE)
#   test_data <- rbind(test_data, sampled_data)
# }
# for (i in seq(from = 0, to = 0.9, by = 0.1)) {
#   lower <- train_pool %>% dplyr::filter(CY3rawperc > YFPrawperc)
#   sampled_data <- lower %>% dplyr::filter(distance > quantile(lower$distance, probs = i) & distance <= quantile(lower$distance, probs = i + 0.1)) %>%
#     dplyr::sample_n(train, replace = FALSE)
#   train_data <- rbind(train_data, sampled_data)
#   lower <- test_pool %>% dplyr::filter(CY3rawperc > YFPrawperc)
#   sampled_data <- lower %>% dplyr::filter(distance > quantile(lower$distance, probs = i) & distance <= quantile(lower$distance, probs = i + 0.1)) %>%
#     dplyr::sample_n(test, replace = FALSE)
#   test_data <- rbind(test_data, sampled_data)
# }


### sampling proportional to quantile range
# radial
train_data = c()
test_data = c()
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  low_q = quantile(all_data$distance, probs = i)[[1]]
  upp_q = quantile(all_data$distance, probs = i + 0.1)[[1]]
  range_data <- all_data %>% dplyr::filter(distance > low_q & distance <= upp_q)
  train = as.integer(nrow(range_data) * (upp_q - low_q))
  sampled_data <- range_data %>% dplyr::sample_n(train + 100, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):nrow(sampled_data),])
}
if (nrow(test_data %>% dplyr::filter(CY3rawperc > 0.85)) == 0) {
  sampled_data <- all_data %>% dplyr::filter(CY3rawperc > 0.85) %>% sample_n(2, replace = FALSE)
  test_data <- rbind(test_data, sampled_data)
}
if (nrow(test_data %>% dplyr::filter(YFPrawperc > 0.85)) == 0) {
  sampled_data <- all_data %>% dplyr::filter(CY3rawperc > 0.85) %>% sample_n(2, replace = FALSE)
  test_data <- rbind(test_data, sampled_data)
}

# biradial
train_data = c()
test_data = c()
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  upper <- all_data %>% dplyr::filter(YFPrawperc > CY3rawperc)
  low_q = quantile(upper$distance, probs = i)[[1]]
  upp_q = quantile(upper$distance, probs = i + 0.1)[[1]]
  range_data <- upper %>% dplyr::filter(distance > low_q & distance <= upp_q)
  train = as.integer(nrow(range_data) * (upp_q - low_q))
  sampled_data <- range_data %>% dplyr::sample_n(train + 50, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):nrow(range_data),])
}
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  lower <- all_data %>% dplyr::filter(CY3rawperc > YFPrawperc)
  low_q = quantile(lower$distance, probs = i)[[1]]
  upp_q = quantile(lower$distance, probs = i + 0.1)[[1]]
  range_data <- lower %>% dplyr::filter(distance > low_q & distance <= upp_q)
  train = as.integer(nrow(range_data) * (upp_q - low_q))
  sampled_data <- range_data %>% dplyr::sample_n(train + 50, replace = FALSE)
  train_data <- rbind(train_data, sampled_data[1:train,])
  test_data <- rbind(test_data, sampled_data[(train + 1):nrow(range_data),])
}

## radial but separate test and train images, log
train_data = c()
test_data = c()
train = 2000
test = 200
for (i in seq(from = 0, to = 0.9, by = 0.1)) {
  sampled_data <- train_pool %>% dplyr::filter(log2distance > quantile(train_pool$log2distance, probs = i) & log2distance <= quantile(train_pool$log2distance, probs = i + 0.1)) %>%
    dplyr::sample_n(train, replace = FALSE)
  train_data <- rbind(train_data, sampled_data)
  sampled_data <- test_pool %>% dplyr::filter(log2distance > quantile(test_pool$log2distance, probs = i) & log2distance <= quantile(test_pool$log2distance, probs = i + 0.1)) %>%
    dplyr::sample_n(test, replace = FALSE)
  test_data <- rbind(test_data, sampled_data)
}



plots = list()
plots[[1]] <- ggplot(train_data, aes(CY3rawperc,YFPrawperc)) + geom_bin2d() + labs(title = "Train Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plots[[2]] <- ggplot(test_data, aes(CY3rawperc,YFPrawperc)) + geom_bin2d() + labs(title = "Test Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plot_grid(plotlist = plots[1:2])

plots = list()
plots[[1]] <- ggplot(train_data, aes(CY3log2trans,YFPlog2trans)) + geom_bin2d() + labs(title = "Train Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plots[[2]] <- ggplot(test_data, aes(CY3log2trans,YFPlog2trans)) + geom_bin2d() + labs(title = "Test Data") + coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
plot_grid(plotlist = plots[1:2])

# use all data
train_data <- test_pool
test_data <- train_pool

# write out the data
write.table(train_data %>% dplyr::select(-contains("CY3"), -contains("CY5"), -contains("YFP")), file = "train_0601_log2_tr3te4.csv", row.names = F, col.names = F, sep = ",")
write.table(train_data %>% dplyr::select(CY3log2trans, YFPlog2trans), file = "train_id_0601_log2_tr3te4.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(-contains("CY3"), -contains("CY5"), -contains("YFP")), file = "test_0601_log2_tr3te4.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(CY3log2trans, YFPlog2trans), file = "test_id_0601_log2_tr3te4.csv", row.names = F, col.names = F, sep = ",")
