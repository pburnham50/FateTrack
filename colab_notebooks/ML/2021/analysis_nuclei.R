# analysis_nuclei.R

# used to look at feature space of nuclei and attribute these to particular markers.

#### Libraries
library(ggplot2); library(umap);
library(tibble); library(dplyr); library(tidyr); library(RANN) ;

setwd("~/Desktop/ml_data/improps_0315/") # set to directory with CSV files containing imaging info

### Functions
getFrame = function(sample, local.radius=100, distal.radius = 500){
  a594.df = data.frame(read.csv(paste0("a594_",sample,".intensity.csv"))) # marker data for a594
  cy3.df = data.frame(read.csv(paste0("cy3_",sample,".intensity.csv"))) # marker data for Cy3
  cy5.df = data.frame(read.csv(paste0("cy5_",sample,".intensity.csv"))) # marker data for Cy5
  gfp.df = data.frame(read.csv(paste0("gfp_",sample,".intensity.csv"))) # marker data for GFP
  
  nuclei = data.frame(read.csv(paste0("nuclei_", sample,".feature.csv"))) # load in nulcear morphology data
  
  merge.df = merge(merge(merge(a594.df, cy3.df, by = "label"), cy5.df, by = "label"), gfp.df, by = "label") # merge all channel data

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
sampleIDS = gsub(pattern = ".intensity.csv", replacement = "",x = gsub(pattern = "gfp_", replacement = "", x = list.files(path = "./",pattern = "gfp.*.csv")))

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
write.csv(subdf, file = "features.csv")

meanInt <- final %>% dplyr::select('meanInt_a594', 'meanInt_cy3', 'meanInt_cy5', 'meanInt_gfp')
write.csv(meanInt, file = "meanInt.csv")

# IGNORE
#subdf_final = sample_n(data.frame(final), 2500)
#subdf =subset.data.frame(subdf_final, select = goodCols)

# scale data and calculate PCA
sca = scale(data.matrix(subdf))
pca <- prcomp(data.matrix(subdf), scale. = TRUE)
plot(pca$x) # plot pca

# calculate umap
umap = umap(sca, n_neighbors = 250, min_dist=0.5)
plot(umap$layout)
umap.df = data.frame(umap$layout)

#alo = cbind(subdf, umap.df)
#mean(log10(alo$meanInt_gfp))

alo = subdf

# find top 5% of gene expression in each channel to annotate cells
alo$secret = (final$meanInt_gfp  > quantile(final$meanInt_gfp,.95))
alo$absorb = (final$meanInt_cy3 > quantile(final$meanInt_cy3,.95))
alo$ta = (final$meanInt_a594 > quantile(final$meanInt_a594,.95))
alo$stem = (final$meanInt_cy5 > quantile(final$meanInt_cy5,.95))

cell_types <- alo %>% dplyr::select("secret", "absorb", "ta", "stem")
write.csv(cell_types, file = "celltypes.csv", row.names = F)

ggplot(alo, aes(X1,X2))+#geom_point(col="lightgrey")+
  geom_point(data = alo[alo$ta,], aes(X1,X2), col = "orange",size=1)+
  geom_point(data = alo[alo$secret,], aes(X1,X2), col = "green",size=1)+
  geom_point(data = alo[alo$absorb,], aes(X1,X2), col = "blue",size=1)+
  geom_point(data = alo[alo$stem,], aes(X1,X2), col = "red",size=1)+
  #scale_color_gradient(high="orange",low="grey")+
  theme_bw()+
  theme(legend.position = "none")

# axis length histograms
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

ggplot(subdf, aes(x = major_axis_length)) + geom_histogram(bins = 100) +
  geom_vline(aes(xintercept=getmode(subdf$major_axis_length)),
             color="blue", linetype="dashed", size=1)

ggplot(subdf, aes(x = minor_axis_length)) + geom_histogram(bins = 100) +
  geom_vline(aes(xintercept=getmode(subdf$minor_axis_length)),
             color="blue", linetype="dashed", size=1)

