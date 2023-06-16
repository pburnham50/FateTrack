library(tibble); library(dplyr); library(tidyr); library(RANN) ; library(tidyverse) ; library(ggplot2); library(cowplot);

setwd("~/Desktop/") # set to directory with CSV files containing imaging info

oldscript <- data.frame(read.csv("sub_0_1_nuclei.feature.csv"))
newscript <- data.frame(read.csv("sub_0_1_nuclei.feature.test.csv"))

colnames(newscript) <- colnames(oldscript)

cols_to_compare <- c("near_contrast", "near_dissim","near_corr", "near_energy", "near_homog", "mid_contrast", "mid_dissim", "mid_corr", "mid_energy", "mid_homog", "far_contrast","far_dissim", "far_corr", "far_energy" 
                                   ,"puncta_count", "puncta_area",  "imvar")

oldscript <- oldscript %>% dplyr::select(cols_to_compare)
newscript <- newscript %>% dplyr::select(cols_to_compare)

oldscript_m <- as.matrix(oldscript)
newscript_m <- as.matrix(newscript)

comp <- oldscript_m - newscript_m

