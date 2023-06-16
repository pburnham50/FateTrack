# Import various libraries
import numpy as np
import time, os, sys
import matplotlib.pyplot as plt
import glob
from scipy import ndimage as ndi
from skimage import color, feature, filters, io, measure, morphology, restoration, segmentation, exposure, restoration, data
import warnings
import math
import pandas as pd
from skimage.feature import greycomatrix, greycoprops
from numpy import genfromtxt

# function
def writeChannelFeatures(imageFrame, channelName, ringMask):
     channel_props = measure.regionprops_table(ringMask, intensity_image=imageFrame, properties=('label','mean_intensity','min_intensity','max_intensity'))

     im_df_2 = pd.DataFrame(channel_props)
     im_df_2 = im_df_2.rename(columns={"mean_intensity": "meanInt_"+channelName, "min_intensity": "minInt_"+channelName, "max_intensity": "maxInt_"+channelName})
     im_df_2["median_intensity_"+channelName] = np.median(imageFrame)
     export = im_df_2.to_csv(output_path + channelName + "_" + sample + ".intensity.csv", index = None, header=True)
     return()


#put your path here
input_path = sys.argv[1] #'~/Desktop/ml_data/data_0527_13dil/'
input_path_seg = sys.argv[2] #'~/Desktop/Trident-main/data_0524_new/'
output_path = sys.argv[3] #'~/Desktop/Trident-main/data_0524_new/'

sample_list = list()
for root, dirs, files in os.walk(input_path, topdown=False):
     for name in dirs:
          sample_list.append(os.path.basename(os.path.normpath(name)))

#sample = sys.argv[4] #
for sample in sample_list:
     newnp = np.asarray(np.load(input_path_seg + sample + '/dapi001_seg.npy',allow_pickle=True)).item()
     outs = (newnp['outlines']>0).astype(int)
     masks = newnp['masks']
     imorig = newnp['img']

     dapi = io.imread( input_path + sample + '/dapi001.tif');
     yfp = io.imread( input_path + sample + '/yfp001.tif');
     #a594 = io.imread( input_path + 'a594_' + sample + '.tif');
     cy3 = io.imread( input_path + sample + '/cy3001.tif');
     cy5 = io.imread( input_path + sample + '/cy5001.tif');



     # initial column names for final output
     column_names = ["near_contrast","near_dissim", "near_corr", "near_energy", "near_homog",
                     "mid_contrast","mid_dissim", "mid_corr", "mid_energy", "mid_homog",
                     "far_contrast","far_dissim", "far_corr", "far_energy", "far_homog",
                     "puncta_count", "puncta_area"]

     image_props = measure.regionprops_table(masks, intensity_image=imorig,properties=('label','area','filled_area',
                                                                                         'bbox', 'centroid',
                                                                                          'eccentricity','solidity','convex_area',
                                                                                          'mean_intensity','min_intensity','max_intensity',
                                                                                          'orientation','major_axis_length','minor_axis_length',
                                                                                          'perimeter','extent','intensity_image'))
     im_df = pd.DataFrame(image_props)
     otherfets = pd.DataFrame()

     contrast_vec = np.empty((1,1))


     for j in range(im_df.shape[0]):
             infoc = np.sum((im_df['intensity_image'][j]-np.mean(im_df['intensity_image'][j]))**2)
             contrast_vec = np.append(contrast_vec,infoc)

             imtest = im_df['intensity_image'][j].astype('uint8')
             minaxis = im_df['minor_axis_length'][j]

             gclm_farrange = greycomatrix(imtest, [minaxis/2], angles=[0], levels=256,
                                     symmetric=False, normed=True)

             gclm_midrange = greycomatrix(imtest, [minaxis/4], angles=[0], levels=256,
                                     symmetric=False, normed=True)

             gclm_closerange = greycomatrix(imtest, distances=[2], angles=[0], levels=256,
                                     symmetric=False, normed=True)

             intimage = im_df['intensity_image'][j]
             spotlabels = measure.label(intimage, connectivity=2, background=0)

             npspots = np.max(np.unique(spotlabels))
             spotpixels = len(spotlabels[spotlabels!=0])


             tmpvals = [greycoprops(gclm_closerange, 'contrast')[0,0],
              greycoprops(gclm_closerange, 'dissimilarity')[0,0],
              greycoprops(gclm_closerange, 'correlation')[0,0],
              greycoprops(gclm_closerange, 'energy')[0,0],
              greycoprops(gclm_closerange, 'homogeneity')[0,0],
              greycoprops(gclm_midrange, 'contrast')[0,0],
              greycoprops(gclm_midrange, 'dissimilarity')[0,0],
              greycoprops(gclm_midrange, 'correlation')[0,0],
              greycoprops(gclm_midrange, 'energy')[0,0],
              greycoprops(gclm_midrange, 'homogeneity')[0,0],
              greycoprops(gclm_farrange, 'contrast')[0,0],
              greycoprops(gclm_farrange, 'dissimilarity')[0,0],
              greycoprops(gclm_farrange, 'correlation')[0,0],
              greycoprops(gclm_farrange, 'energy')[0,0],
              greycoprops(gclm_farrange, 'homogeneity')[0,0],
              npspots, spotpixels]

             erdf = pd.DataFrame(tmpvals).T

             otherfets = otherfets.append(erdf, ignore_index=True)

     otherfets.columns = column_names
     #    df_c = pd.concat([im_df.reset_index(drop=True), xs], axis=1)
     df_c = pd.concat([im_df.reset_index(drop=True), otherfets], axis=1)
     tmpdf = pd.DataFrame(contrast_vec[1:])
     tmpdf.columns = ['imvar']

     im2_df = pd.concat([df_c.reset_index(drop=True), tmpdf], axis=1)


     ### look in anulus around nuclei to detect channel intensity
     # change these to change the ring of pixels around the nucleus
     structmax = morphology.square(11)
     structmin = morphology.square(3)

     ringmask = morphology.dilation(masks,structmax)-morphology.dilation(masks,structmin)

     #writeChannelFeatures(yfp, "yfp", ringMask = ringmask)
     #writeChannelFeatures(a594, "a594", ringMask = ringmask)
     #writeChannelFeatures(cy3, "cy3", ringMask = ringmask)
     #writeChannelFeatures(cy5, "cy5", ringMask = ringmask)

     export = im2_df.to_csv(output_path + sample + '/' + sample + "_nuclei.feature.csv", index = None, header=True)
     #export = im_df_2.to_csv(output_path + sample + ".red.intensity.csv", index = None, header=True)
