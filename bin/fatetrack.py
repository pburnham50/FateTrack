from __future__ import print_function
import numpy as np
import time, os, sys
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage import color, feature, filters, io, measure, morphology, segmentation, img_as_ubyte, transform
import warnings
import math
import pandas as pd
import argparse
import subprocess
import re
import glob
from skimage.segmentation import clear_border
from ortools.graph import pywrapgraph
import time

def buildFeatureFrame(filename,timepoint,pathtoimage="./"):
    temp = np.asarray(np.load(filename,allow_pickle=True)).item()
    imfilename = temp['filename'].split('/')[-1]
    img = io.imread(pathtoimage+imfilename);
    masks = clear_border(temp['masks'])
    image_props = measure.regionprops_table(masks,
                                        intensity_image=img,
                                        properties=('label','area','filled_area', 'centroid',
                                                    'eccentricity','mean_intensity'))
    im_df = pd.DataFrame(image_props)
    im_df['time'] = timepoint
    return(im_df)


def buildOffsetFrame(filename_t0,filename_t1,pathtoimage="./"):

    temp1 = np.asarray(np.load(filename_t0,allow_pickle=True)).item()
    imfilename1 = temp1['filename'].split('/')[-1]
    img1 = io.imread(pathtoimage+imfilename1);

    temp2 = np.asarray(np.load(filename_t1,allow_pickle=True)).item()
    imfilename2 = temp2['filename'].split('/')[-1]
    img2 = io.imread(pathtoimage+imfilename2);

    masks = clear_border(temp1['masks'])

    image_props = measure.regionprops_table(temp1['masks'],
                                        intensity_image=img2,
                                        properties=('label','area','filled_area', 'centroid',
                                                    'eccentricity','mean_intensity'))
    im_df = pd.DataFrame(image_props)
    im_df['time'] = None
    return(im_df)
