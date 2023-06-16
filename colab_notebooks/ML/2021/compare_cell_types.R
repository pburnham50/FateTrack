# script to compare cell types identified by analysis_nuclei.R when size of nuclear ring is modified in feature_ext.HCR.py

library(dplyr);

setwd("~/Desktop/data/data1008/improps/")

### comparison of mean intensity
# load in files
max11_min3 <- read.csv("max11min3/meanInt_max11min3.csv")
max11_min5 <- read.csv("max11min5/meanInt_max11min5.csv")
max16_min3 <- read.csv("max16min3/meanInt_max16min3.csv")
max16_min5 <- read.csv("max16min5/meanInt_max16min5.csv")

# create comparison tables with max11_min3 as baseline
comparison_df <- function (parameters) {
  new_df <- ((parameters$'meanInt_a594' - max11_min3$'meanInt_a594') / (max11_min3$'meanInt_a594')) %>% as.data.frame()
  for (stain in c('meanInt_cy3', 'meanInt_cy5', 'meanInt_gfp')) {
    new_df <- new_df %>% cbind((parameters$stain - max11_min3$stain) / (max11_min3$stain))
  }
  new_df
}
compare_max11min5 <- comparison_df(max11_min5)


### comparison of which cells are in top 5% of expression
# load in files
max11_min3 <- read.csv("max11min3/celltypes_max11min3.csv")
max11_min5 <- read.csv("max11min5/celltypes_max11min5.csv")
max16_min3 <- read.csv("max16min3/celltypes_max16min3.csv")
max16_min5 <- read.csv("max16min5/celltypes_max16min5.csv")

# combine by cell type
secret <- cbind(max11_min3$secret, max11_min5$secret, max16_min3$secret, max16_min5$secret) %>% as.data.frame()
absorb <- cbind(max11_min3$absorb, max11_min5$absorb, max16_min3$absorb, max16_min5$absorb) %>% as.data.frame()
ta <- cbind(max11_min3$ta, max11_min5$ta, max16_min3$ta, max16_min5$ta) %>% as.data.frame()
stem <- cbind(max11_min3$stem, max11_min5$stem, max16_min3$stem, max16_min5$stem) %>% as.data.frame()

# rename columns by parameters used
new_colnames <- c('max11min3', 'max11min5', 'max16min3', 'max16min5')
names(secret) <- new_colnames
names(absorb) <- new_colnames
names(ta) <- new_colnames
names(stem) <- new_colnames

# calculate matches between other parameters and baseline of max11min3
secret <- secret %>% dplyr::mutate(compare_max11min5 = secret$max11min3 == secret$max11min5,
                                   compare_max16min3 = secret$max11min3 == secret$max16min3,
                                   compare_max16min5 = secret$max11min3 == secret$max16min5)
absorb <- absorb %>% dplyr::mutate(compare_max11min5 = absorb$max11min3 == absorb$max11min5,
                                   compare_max16min3 = absorb$max11min3 == absorb$max16min3,
                                   compare_max16min5 = absorb$max11min3 == absorb$max16min5)
ta <- ta %>% dplyr::mutate(compare_max11min5 = ta$max11min3 == ta$max11min5,
                                   compare_max16min3 = ta$max11min3 == ta$max16min3,
                                   compare_max16min5 = ta$max11min3 == ta$max16min5)
stem <- stem %>% dplyr::mutate(compare_max11min5 = stem$max11min3 == stem$max11min5,
                                   compare_max16min3 = stem$max11min3 == stem$max16min3,
                                   compare_max16min5 = stem$max11min3 == stem$max16min5)

# calculate difference
difference <- function(celltype, compare_parameter) {
  c <- 0
  for (r in celltype[,compare_parameter]) {
    if (!r) {
      c = c + 1
    }
  }
  c
}

# create data frame for each nucleus settings
compare_max11min5 <- c(difference(secret, 'compare_max11min5'), difference(absorb, 'compare_max11min5'), difference(ta, 'compare_max11min5'), difference(stem, 'compare_max11min5')) 
dim(compare_max11min5) <- c(1,4)
compare_max11min5 <- compare_max11min5 %>% as.data.frame()
compare_max16min3 <- c(difference(secret, 'compare_max16min3'), difference(absorb, 'compare_max16min3'), difference(ta, 'compare_max16min3'), difference(stem, 'compare_max16min3'))
dim(compare_max16min3) <- c(1,4)
compare_max16min3 <- compare_max16min3 %>% as.data.frame()
compare_max16min5 <- c(difference(secret, 'compare_max16min5'), difference(absorb, 'compare_max16min5'), difference(ta, 'compare_max16min5'), difference(stem, 'compare_max16min5'))
dim(compare_max16min5) <- c(1,4)
compare_max16min5 <- compare_max16min5 %>% as.data.frame()

# change column names
new_colnames <- c('secret', 'absorb', 'ta', 'stem')
names(compare_max11min5) <- new_colnames
names(compare_max16min3) <- new_colnames
names(compare_max16min5) <- new_colnames

