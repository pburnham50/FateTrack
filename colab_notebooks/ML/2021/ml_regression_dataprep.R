# adapted from ml_celltype_dataprep.R

# used to look at feature space of nuclei and attribute these to particular celltype markers.
# prepares the data for use in machine learning.
# specifically prepares data for regression on features
# the difference is in how each cell is assigned a "percent" of a certain celltype

#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ;

setwd("~/Desktop/ml_data/data_0415/") # set to directory with CSV files containing imaging info

### Functions
getFrame = function(sample, local.radius=100, distal.radius = 500){
  a594.df = data.frame(read.csv(paste0("a594_",sample,".intensity.csv"))) # marker data for a594
  cy3.df = data.frame(read.csv(paste0("cy3_",sample,".intensity.csv"))) # marker data for Cy3
  cy5.df = data.frame(read.csv(paste0("cy5_",sample,".intensity.csv"))) # marker data for Cy5
  yfp.df = data.frame(read.csv(paste0("yfp_",sample,".intensity.csv"))) # marker data for GFP
  
  nuclei = data.frame(read.csv(paste0("nuclei_", sample,".feature.csv"))) # load in nulcear morphology data
  
  merge.df = merge(merge(merge(a594.df, cy3.df, by = "label"), cy5.df, by = "label"), yfp.df, by = "label") # merge all channel data
  
  merge.final = merge(nuclei, merge.df, by = "label")
  merge.final = merge.final[!(colnames(merge.final) %in% c("intensity.image")),]
  merge.final$sample = sample
  merge.final = merge.final[complete.cases(merge.final), ]
  
  merge.final$ratio = abs(merge.final$bbox.0 - merge.final$bbox.2)/abs(merge.final$bbox.1 - merge.final$bbox.3) # calculate aspect ratio
  # calculate minimum distance (centroid to centroid) between nuclei. also number of cells in a given radius.
  sub.tab = data.frame(subset(merge.final, select = c("centroid.0","centroid.1")))
  if (nrow(merge.final)<2){
    merge.final$mindist = 10**5
    merge.final$nearby.cells = 0
    merge.final$mid.cells = 0
    
  }else{
    # look at distances between every pair of nuclei and take the minimum
    mindist = apply(data.matrix(as.matrix(dist(sub.tab))),1,FUN=function(x) {min(x[x > 0])})
    merge.final$mindist = mindist
    
    # treat every nucleus as a point (centroid), count all other nuclei within given radius
    nnframe = nn2(data = sub.tab,radius = local.radius, searchtype = 'radius', k = nrow(sub.tab))
    merge.final$nearby.cells = rowSums(data.matrix(nnframe$nn.idx[,-1]>0))
    
    nnframe = nn2(data = sub.tab,radius = local.radius, searchtype = 'radius', k = nrow(sub.tab))
    merge.final$mid.cells = rowSums(data.matrix(nnframe$nn.idx[,-1]>0))
  }
  
  #dtab$mean_intensity = dtab$mean_intensity/median(dtab$mean_intensity)
  
  
  return(merge.final)
}

### Execution

# get all image data names
sampleIDS = gsub(pattern = ".intensity.csv", replacement = "",x = gsub(pattern = "yfp_", replacement = "", x = list.files(path = "./",pattern = "yfp.*.csv")))

# construct large dataframe containing data for all images
final = c()

for (j in sampleIDS){
  final = rbind(final, getFrame(j))
}

# filter out very large nuclei
final = (final[complete.cases(final), ])
final = final[final$area < 3000,]
final = final[final$mindist > 0,]

# these are the columns we use for PCA/UMAP
goodCols = c("area", "filled_area","eccentricity", "solidity", "convex_area" #"mean_intensity", "min_intensity", "max_intensity",
             ,"orientation", "major_axis_length","minor_axis_length", "perimeter", "extent" 
             ,"near_contrast", "near_dissim","near_corr", "near_energy", "near_homog", "mid_contrast", "mid_dissim", "mid_corr", "mid_energy", "mid_homog", "far_contrast","far_dissim", "far_corr", "far_energy" 
             ,"puncta_count", "puncta_area", "mindist", "nearby.cells", "mid.cells", "ratio")

subdf<-subset.data.frame(final, select = goodCols)
#write.csv(subdf, file = "features.csv")

meanInt <- final %>% dplyr::select('meanInt_a594', 'meanInt_cy3', 'meanInt_cy5', 'meanInt_yfp')
#write.csv(meanInt, file = "meanInt.csv")


feat_type = subdf


# original normalization
feat_type$CY3 = (final$meanInt_cy3)#/min(feat_type$CY3_background))
feat_type$CY5 = (final$meanInt_cy5)#/min(feat_type$CY5_background))
feat_type$YFP = (final$meanInt_yfp)#/min(feat_type$YFP_background))
feat_type$CY3minnorm = (feat_type$CY3/min(feat_type$CY3))
feat_type$CY5minnorm = (feat_type$CY5/min(feat_type$CY5))
feat_type$YFPminnorm = (feat_type$YFP/min(feat_type$YFP))
feat_type$CY3norm = feat_type$CY3minnorm / (feat_type$CY3minnorm + feat_type$CY5minnorm + feat_type$YFPminnorm)
feat_type$CY5norm = feat_type$CY5minnorm / (feat_type$CY3minnorm + feat_type$CY5minnorm + feat_type$YFPminnorm)
feat_type$YFPnorm = feat_type$YFPminnorm / (feat_type$CY3minnorm + feat_type$CY5minnorm + feat_type$YFPminnorm)

#plot(feat_type$CY5norm, feat_type$YFPnorm)

# re-normalization by percent
feat_type$CY3perc = ((feat_type$CY3norm - min(feat_type$CY3norm))/(max(feat_type$CY3norm) - min(feat_type$CY3norm)))
feat_type$CY5perc = ((feat_type$CY5norm - min(feat_type$CY5norm))/(max(feat_type$CY5norm) - min(feat_type$CY5norm)))
feat_type$YFPperc = ((feat_type$YFPnorm - min(feat_type$YFPnorm))/(max(feat_type$YFPnorm) - min(feat_type$YFPnorm)))
feat_type$CY3normperc = (feat_type$CY3perc / (feat_type$CY3perc + feat_type$CY5perc + feat_type$YFPperc))
feat_type$CY5normperc = (feat_type$CY5perc / (feat_type$CY3perc + feat_type$CY5perc + feat_type$YFPperc))
feat_type$YFPnormperc = (feat_type$YFPperc / (feat_type$CY3perc + feat_type$CY5perc + feat_type$YFPperc))

#plot(feat_type$CY5normperc, feat_type$YFPnormperc)


# normalization on raw intensities
feat_type$CY3rawperc = ((feat_type$CY3 - min(feat_type$CY3))/(max(feat_type$CY3) - min(feat_type$CY3)))
feat_type$CY5rawperc = ((feat_type$CY5 - min(feat_type$CY5))/(max(feat_type$CY5) - min(feat_type$CY5)))
feat_type$YFPrawperc = ((feat_type$YFP - min(feat_type$YFP))/(max(feat_type$YFP) - min(feat_type$YFP)))
feat_type$CY3rawnormperc = (feat_type$CY3rawperc / (feat_type$CY3rawperc + feat_type$CY5rawperc + feat_type$YFPrawperc))
feat_type$CY5rawnormperc = (feat_type$CY5rawperc / (feat_type$CY3rawperc + feat_type$CY5rawperc + feat_type$YFPrawperc))
feat_type$YFPrawnormperc = (feat_type$YFPrawperc / (feat_type$CY3rawperc + feat_type$CY5rawperc + feat_type$YFPrawperc))

#plot(feat_type$YFPnormperc, feat_type$CY5normperc)
hist(feat_type$CY3perc)

# prepare data for ml
sampled_data <- feat_type %>% dplyr::sample_n(2200, replace = FALSE)
train_data <- sampled_data[1:2000,]
test_data <- sampled_data[2001:2200,]


write.table(train_data %>% dplyr::select(-contains("CY3"), -contains("CY5"), -contains("YFP")), file = "train_0520_raw.csv", row.names = F, col.names = F, sep = ",")
write.table(train_data %>% dplyr::select(CY5rawnormperc, YFPrawnormperc), file = "train_id_0520_raw.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(-contains("CY3"), -contains("CY5"), -contains("YFP")), file = "test_0520_raw.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(CY5rawnormperc, YFPrawnormperc), file = "test_id_0520_raw.csv", row.names = F, col.names = F, sep = ",")
