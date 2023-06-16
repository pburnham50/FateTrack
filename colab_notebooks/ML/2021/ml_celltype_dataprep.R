# adapted from analysis_nuclei.R

# used to look at feature space of nuclei and attribute these to particular celltype markers.
# prepares the data for use in machine learning.

#### Libraries
library(tibble); library(dplyr); library(tidyr); library(RANN) ;

setwd("~/Desktop/ml_data/HCR_0610/507_D3_1_w1/") # set to directory with CSV files containing imaging info

### Functions
getFrame = function(sample, local.radius=100, distal.radius = 500){
  #a594.df = data.frame(read.csv(paste0("a594_",sample,".intensity.csv"))) # marker data for a594
  #cy3.df = data.frame(read.csv(paste0("cy3_",sample,".intensity.csv"))) # marker data for Cy3
  #cy5.df = data.frame(read.csv(paste0("cy5_",sample,".intensity.csv"))) # marker data for Cy5
  #yfp.df = data.frame(read.csv(paste0("yfp_",sample,".intensity.csv"))) # marker data for yfp
  
  #nuclei = data.frame(read.csv(paste0("nuclei_", sample,".feature.csv"))) # load in nulcear morphology data
  
  cy3.df = data.frame(read.csv(paste0(sample,"2CY3_measurements.csv"))) # marker data for Cy3
  cy5.df = data.frame(read.csv(paste0(sample,"2CY5_measurements.csv"))) # marker data for Cy5
  yfp.df = data.frame(read.csv(paste0(sample,"2YFP_measurements.csv"))) # marker data for GFP
  
  nuclei = data.frame(read.csv(paste0(sample,"2_nuclei.feature.csv"))) # load in nulcear morphology data
  
  
  merge.df = merge(merge(merge(cy3.df, cy5.df, by = c("subimage", "label")), yfp.df, by = c("subimage", "label"))) # merge all channel data
  
  merge.final = merge(nuclei, merge.df, by = c("subimage", "label")
  merge.final = merge.final[!(colnames(merge.final) %in% c("intensity.image")),]
  merge.final$sample = sample
  merge.final = merge.final[complete.cases(merge.final), ]
  
  merge.final$ratio = abs(merge.final$bbox.0 - merge.final$bbox.2)/abs(merge.final$bbox.1 - merge.final$bbox.3) # calculate aspect ratio
  # calculate minimum distance (centroid to centroid) between nuclei. also number of cells in a given radius.
  # sub.tab = data.frame(subset(merge.final, select = c("centroid.0","centroid.1")))
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
  # 
  #dtab$mean_intensity = dtab$mean_intensity/median(dtab$mean_intensity)
  
  
  return(merge.final)
}

### Execution

# get all image data names
#sampleIDS = gsub(pattern = ".intensity.csv", replacement = "",x = gsub(pattern = "yfp_", replacement = "", x = list.files(path = "./",pattern = "yfp.*.csv")))
sampleIDS = gsub(pattern = "_measurements.csv", replacement = "",x = gsub(pattern = "_yfp", replacement = "", x = list.files(path = "./",pattern = "*_yfp_measurements.csv")))

# construct large dataframe containing data for all images
final = c()

for (j in sampleIDS){
  print(j)
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

meanInt <- final %>% dplyr::select('meanInt_cy3', 'meanInt_cy5', 'meanInt_yfp')
#write.csv(meanInt, file = "meanInt.csv")


feat_type = subdf

# THRESHOLD APPROACH find top 5% of gene expression in each channel to annotate cells
# feat_type$enterocyte = (final$meanInt_yfp  > quantile(final$meanInt_yfp,.99))
# feat_type$paneth = (final$meanInt_cy3 > quantile(final$meanInt_cy3,.99))
# feat_type$ta = (final$meanInt_a594 > quantile(final$meanInt_a594,.99))
# feat_type$goblet = (final$meanInt_cy5 > quantile(final$meanInt_cy5,.99))

#cell_types <- feat_type %>% dplyr::select("secret", "absorb", "ta", "stem")
#write.csv(cell_types, file = "celltypes.csv", row.names = F)

#write.csv(feat_type, file = "features_with_celltypes.csv", row.names = F)

# TOP N NUCLEI APPROACH: find top 100 of each cell type
# yfp <- sort(final$meanInt_yfp, decreasing = TRUE)
# yfp_n <- as.integer(yfp[100])
# feat_type$enterocyte = (final$meanInt_yfp > yfp_n)
# cy5 <- sort(final$meanInt_cy5, decreasing = TRUE)
# cy5_n <- as.integer(cy5[100])
# feat_type$goblet = (final$meanInt_cy5 > cy5_n)

# WEIGHTED EXPRESSION APPROACH
feat_type$CY3 = (final$meanInt_cy3)#/min(feat_type$CY3_background))
feat_type$CY5 = (final$meanInt_cy5)#/min(feat_type$CY5_background))
feat_type$YFP = (final$meanInt_yfp)#/min(feat_type$YFP_background))
feat_type$CY3 = (feat_type$CY3/min(feat_type$CY3))
feat_type$CY5 = (feat_type$CY5/min(feat_type$CY5))
feat_type$YFP = (feat_type$YFP/min(feat_type$YFP))
feat_type$CY3norm = feat_type$CY3 / (feat_type$CY3 + feat_type$CY5 + feat_type$YFP)
feat_type$CY5norm = feat_type$CY5 / (feat_type$CY3 + feat_type$CY5 + feat_type$YFP)
feat_type$YFPnorm = feat_type$YFP / (feat_type$CY3 + feat_type$CY5 + feat_type$YFP)
feat_type$magnitude =  sqrt((feat_type$YFP)**2 + (feat_type$CY5)**2 + (feat_type$CY3)**2)
feat_type_high <- feat_type %>% dplyr::filter(feat_type$magnitude>quantile(feat_type$magnitude,probs = .85))
#feat_type_low = feat_type[feat_type$magnitude<quantile(feat_type$magnitude,probs = .15),]
feat_type_high <- feat_type_high %>% dplyr::mutate(enterocyte = YFPnorm>quantile(YFPnorm,probs = .85),
                                                   clu = CY3norm > quantile(CY3norm, probs = 0.85))
feat_type <- feat_type_high
plot(feat_type$YFPnorm, feat_type$CY3norm)
plot(feat_type_high$YFPnorm, feat_type_high$CY3norm)


# only want goblet and enterocytes
#feat_type_filtered <- feat_type %>% dplyr::filter(enterocyte == TRUE | goblet == TRUE)
feat_type_filtered <- feat_type %>% dplyr::filter(xor(enterocyte == TRUE, clu == TRUE))
# type column: enterocyte = 1, goblet = 0
feat_type_filtered <- feat_type_filtered %>% dplyr::mutate(type = as.numeric(feat_type_filtered$enterocyte))
enterocytes <- feat_type_filtered %>% dplyr::filter(type == 1)
clu <- feat_type_filtered %>% dplyr::filter(type == 0)



enterocytes_sample <- enterocytes %>% dplyr::sample_n(250, replace = FALSE)
clu_sample <- clu %>% dplyr::sample_n(250, replace = FALSE)

# prepare data for ml
train_data <- dplyr::bind_rows(enterocytes_sample[1:200,], clu_sample[1:200,])
test_data <- dplyr::bind_rows(enterocytes_sample[201:250,], clu_sample[201:250,])

# #for some reason it won't let me remove the column names?
# #write.table(train_data %>% dplyr::select(-enterocyte, - paneth, - ta, -goblet, -type), file = "train_0401.csv", row.names = F, col.names = F, sep = ",")
# #write.table(train_data %>% dplyr::select(type), file = "train_id_0401.csv", row.names = F, col.names = F, sep = ",")
# #write.table(test_data %>% dplyr::select(-enterocyte, - paneth, - ta, -goblet, -type), file = "test_0401.csv", row.names = F, col.names = F, sep = ",")
# #write.table(test_data %>% dplyr::select(type), file = "test_id_0401.csv", row.names = F, col.names = F, sep = ",")
# 
# # write.table(train_data %>% dplyr::select(-enterocyte, -goblet, -type), file = "train_0406.csv", row.names = F, col.names = F, sep = ",")
# # write.table(train_data %>% dplyr::select(type), file = "train_id_0406.csv", row.names = F, col.names = F, sep = ",")
# # write.table(test_data %>% dplyr::select(-enterocyte, -goblet, -type), file = "test_0406.csv", row.names = F, col.names = F, sep = ",")
# # write.table(test_data %>% dplyr::select(type), file = "test_id_0406.csv", row.names = F, col.names = F, sep = ",")
# 
# 
write.table(train_data %>% dplyr::select(-contains("CY3"), -contains("CY5"), -contains("YFP"), -magnitude, -enterocyte, -clu, -type), file = "train_0528_rf.csv", row.names = F, col.names = F, sep = ",")
write.table(train_data %>% dplyr::select(type), file = "train_id_0528_rf.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(-contains("CY3"), -contains("CY5"), -contains("YFP"), -magnitude, -enterocyte, -clu, -type), file = "test_0528_rf.csv", row.names = F, col.names = F, sep = ",")
write.table(test_data %>% dplyr::select(type), file = "test_id_0528_rf.csv", row.names = F, col.names = F, sep = ",")

#just top features
# write.table(train_data %>% dplyr::select(major_axis_length, eccentricity, perimeter, ratio, convex_area,
#                                          near_energy, nearby.cells, near_homog, mid.cells), file = "train_0423.csv", row.names = F, col.names = F, sep = ",")
# write.table(train_data %>% dplyr::select(type), file = "train_id_0423.csv", row.names = F, col.names = F, sep = ",")
# write.table(test_data %>% dplyr::select(major_axis_length, eccentricity, perimeter, ratio, convex_area,
#                                         near_energy, nearby.cells, near_homog, mid.cells), file = "test_0423.csv", row.names = F, col.names = F, sep = ",")
# write.table(test_data %>% dplyr::select(type), file = "test_id_0423.csv", row.names = F, col.names = F, sep = ",")

# Data prep without splitting (python script will split into train and test)
#ml_data <- dplyr::bind_rows(enterocytes_sample, goblet_sample)
#write.table(ml_data %>% dplyr::select(-CY3, -CY5, -YFP, -CY3norm, -CY5norm, -YFPnorm, -magnitude, -enterocyte, -goblet, -type), file = "data_0413.csv", row.names = F, col.names = F, sep = ",")
#write.table(ml_data %>% dplyr::select(type), file = "target_0413.csv", row.names = F, col.names = F, sep = ",")
