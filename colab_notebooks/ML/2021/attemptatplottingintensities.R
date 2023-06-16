library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(ggplot2)

setwd('~/Desktop/ml_data/HCR/HCR_508_D5_1_w1/subimages/')

yfp_intensity <- data.frame(read.csv('sub_0_2CY3_measurements.csv'))
yfp_intensity$'X2CY3_max_intensity'=as.numeric(levels(yfp_intensity$'X2CY3_max_intensity'))[yfp_intensity$'X2CY3_max_intensity']

ggplot() + geom_point(data = yfp_intensity, aes(x = centroid.1, y = -centroid.0, col = 'X2CY3_max_intensity')) +
  scale_color_continuous(low = "purple", high = "yellow")+theme_dark()+theme(legend.position = "none")

cy5_intensity <- data.frame(read.csv('/Users/laila/Desktop/Trident-main/data_0520/403_1/403_1_cy5_meaurements.csv'))
ggplot() + geom_point(data = cy5_intensity, aes(x = centroid.1, y = -centroid.0, col = cy5_mean_intensity)) +
  scale_color_continuous(low = "blue", high = "yellow")+theme_dark()+theme(legend.position = "none")
