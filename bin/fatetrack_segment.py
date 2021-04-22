"""
Philip Burnham
Arjun Raj Lab
FateTrack - registration
2021
fatetrack_segment.py
"""

import matplotlib
import matplotlib.pyplot as plt
import time, os, sys
import pandas as pd
import numpy as np
import glob
import skimage
from skimage import io,registration,transform,filters,img_as_uint,exposure,restoration
import cellpose
from cellpose import utils,models

def segmentCellpose(inputImage, cellpose_out, nucDiameter=0.,cellprob_threshold=-3.5, flow_threshold=0.6):
    #multitiff = io.MultiImage(inputImage)
    #image_series = io.collection.concatenate_images(multitiff)

    image_series = io.imread(inputImage)

    if  1-(os.path.isdir(cellpose_out)):
            os.mkdir(cellpose_out)

    for j in range(image_series.shape[0]):
        skimage.io.imsave(fname=cellpose_out+"/"+str(j)+".tif",arr=image_series[j])


    model = models.Cellpose(gpu=False, model_type='nuclei')
    image_series_list = []
    [image_series_list.append(exposure.equalize_adapthist(image_series[i])) for i in range(image_series.shape[0])]

    masks, flows, styles, diams = model.eval(image_series_list, diameter=nucDiameter, \
                                            cellprob_threshold=cellprob_threshold, flow_threshold=flow_threshold, channels=[0,0])

    files = [cellpose_out+"/"+str(x)+".tif" for x in range(image_series.shape[0])]
    cellpose.io.masks_flows_to_seg(image_series_list, masks, flows, diams, files)


    return(masks)
