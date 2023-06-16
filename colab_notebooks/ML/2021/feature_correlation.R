setwd("~/Desktop/ml_data/HCR_0610/507_D3_1_w1/")

n1 <- data.frame(read.csv("1_nuclei.feature.csv"))
n2 <- data.frame(read.csv("2_nuclei.feature.csv"))

colnames(n1) <- c("subimage", "label", "area", "filled_area",
                      "bbox.0", "bbox.1", "bbox.2", "bbox.3", "centroid.0", "centroid.1",
                      "eccentricity", "solidity", "convex_area", "mean_intensity", "min_intensity", "max_intensity",
                      "orientation", "major_axis_length", "minor_axis_length", "perimeter", "extent", "intensity.image",
                      "near_contrast","near_dissim", "near_corr", "near_energy", "near_homog",
                      "mid_contrast","mid_dissim", "mid_corr", "mid_energy", "mid_homog",
                      "far_contrast","far_dissim", "far_corr", "far_energy", "far_homog",
                      "puncta_count", "puncta_area", "imvar")
colnames(n2) <- c("subimage", "label", "area", "filled_area",
                      "bbox.0", "bbox.1", "bbox.2", "bbox.3", "centroid.0", "centroid.1",
                      "eccentricity", "solidity", "convex_area", "mean_intensity", "min_intensity", "max_intensity",
                      "orientation", "major_axis_length", "minor_axis_length", "perimeter", "extent", "intensity.image",
                      "near_contrast","near_dissim", "near_corr", "near_energy", "near_homog",
                      "mid_contrast","mid_dissim", "mid_corr", "mid_energy", "mid_homog",
                      "far_contrast","far_dissim", "far_corr", "far_energy", "far_homog",
                      "puncta_count", "puncta_area", "imvar")

n1 <- dplyr::semi_join(n1, n2, by = c("subimage", "label"))

n1_m <- data.matrix(n1)
n2_m <- data.matrix(n2)
nf_diff <- matrix(nrow = 115648, ncol = 40)
nf_diff[,1:2] <- n1_m[,1:2]
nf_diff[,5:10] <- n1_m[,5:10]
for (i in c(3, 4, 11:40)) {
  nf_diff[,i] <- n1_m[,i] - n2_m[,i]
}
colnames(nf_diff) <- c("subimage", "label", "area", "filled_area",
                       "bbox.0", "bbox.1", "bbox.2", "bbox.3", "centroid.0", "centroid.1",
                       "eccentricity", "solidity", "convex_area", "mean_intensity", "min_intensity", "max_intensity",
                       "orientation", "major_axis_length", "minor_axis_length", "perimeter", "extent", "intensity.image",
                       "near_contrast","near_dissim", "near_corr", "near_energy", "near_homog",
                       "mid_contrast","mid_dissim", "mid_corr", "mid_energy", "mid_homog",
                       "far_contrast","far_dissim", "far_corr", "far_energy", "far_homog",
                       "puncta_count", "puncta_area", "imvar")

featureCols = c("area", "filled_area","eccentricity", "solidity", "convex_area" #"mean_intensity", "min_intensity", "max_intensity",
                ,"orientation", "major_axis_length","minor_axis_length", "perimeter", "extent" 
                ,"near_contrast", "near_dissim","near_corr", "near_energy", "near_homog", "mid_contrast", "mid_dissim", "mid_corr", "mid_energy", "mid_homog", "far_contrast","far_dissim", "far_corr", "far_energy" 
                ,"puncta_count", "puncta_area")

featureCols = c("centroid.0", "centroid.1")

n1_feat <- n1 %>% dplyr::select(featureCols)
n2_feat <- n2 %>% dplyr::select(featureCols)
cor_same <- cor(n1_feat, n2_feat)

n2_feat_shuffle <- n2_feat %>% sample_frac(1, replace = F)
cor_random <- cor(n1_feat, n2_feat_shuffle)

correlations <- matrix(nrow = 2, ncol = 2)
for(i in seq(1, 2)) {
  correlations[i, 1] <- cor_same[i, i]
  correlations[i, 2] <- cor_random[i, i]
}
rownames(correlations) = featureCols
