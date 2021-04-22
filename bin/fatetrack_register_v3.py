"""
Philip Burnham
Arjun Raj Lab
FateTrack - registration
2021
fatetrack_register.py
"""

import matplotlib
import matplotlib.pyplot as plt
import time, os, sys
import pandas as pd
import numpy as np
import glob
from skimage import io,registration,transform,filters,img_as_uint,exposure,restoration

def register(image_series, additionalCropPxl):
    numFrames = image_series.shape[0]
    shift_vectors = [registration.phase_cross_correlation(image_series[i], image_series[i+1])[0] for i in range((numFrames-1))]
    shift_vectors.insert(0, np.asarray([0.0,0.0]))

    # set null values for later assignment
    image_registered = [np.zeros(image.shape) for image in image_series]
    cumm_shift_row,cumm_shift_col=0,0
    max_cumm_shift_row,max_cumm_shift_col=0,0
    min_cumm_shift_row,min_cumm_shift_col=0,0

    #replace frames with appropriate translations
    for i in range(len(image_series)):
        shift_row = int(shift_vectors[i][0])
        shift_col = int(shift_vectors[i][1])
        cumm_shift_row = cumm_shift_row+shift_row
        cumm_shift_col = cumm_shift_col+shift_col

        max_cumm_shift_row = np.max((max_cumm_shift_row,(cumm_shift_row)))
        max_cumm_shift_col = np.max((max_cumm_shift_col,(cumm_shift_col)))
        min_cumm_shift_row = np.min((min_cumm_shift_row,(cumm_shift_row)))
        min_cumm_shift_col = np.min((min_cumm_shift_col,(cumm_shift_col)))

        tmp_trform = transform.AffineTransform(translation=[-cumm_shift_col,-cumm_shift_row])
        image_registered[i] = transform.warp(image_series[i], tmp_trform)

    # crop images to remove black borders.
    imShape = image_registered[0].shape
    cropList = []
    for j in range(len(image_series)):
        cropped = image_registered[j][(max_cumm_shift_row+additionalCropPxl):(imShape[0]-(np.abs(min_cumm_shift_row)+additionalCropPxl)),(max_cumm_shift_col+additionalCropPxl):(imShape[1]-(np.abs(min_cumm_shift_col)+additionalCropPxl))]
        cropped = img_as_uint(cropped/np.max(cropped))
        cropList.append(cropped)
    cropList = np.stack(cropList,axis = 0)

    shiftVecDF = pd.DataFrame(shift_vectors)
    shiftVecDF.columns = ('shiftRow','shiftColumn')

    return(cropList,shiftVecDF)


# add in definition to take a shift vector and apply transforms to an image.
def applyRegistration(image_series, shiftFrame,additionalCropPxl):
    numFrames = image_series.shape[0]
    shift_vectors = np.asarray(shiftFrame)

    image_registered = [np.zeros(image.shape) for image in image_series]
    cumm_shift_row,cumm_shift_col=0,0
    max_cumm_shift_row,max_cumm_shift_col=0,0
    min_cumm_shift_row,min_cumm_shift_col=0,0

    #replace frames with appropriate translations
    for i in range(len(image_series)):
        shift_row = int(shift_vectors[i][0])
        shift_col = int(shift_vectors[i][1])
        cumm_shift_row = cumm_shift_row+shift_row
        cumm_shift_col = cumm_shift_col+shift_col

        max_cumm_shift_row = np.max((max_cumm_shift_row,(cumm_shift_row)))
        max_cumm_shift_col = np.max((max_cumm_shift_col,(cumm_shift_col)))
        min_cumm_shift_row = np.min((min_cumm_shift_row,(cumm_shift_row)))
        min_cumm_shift_col = np.min((min_cumm_shift_col,(cumm_shift_col)))

        tmp_trform = transform.AffineTransform(translation=[-cumm_shift_col,-cumm_shift_row])
        image_registered[i] = transform.warp(image_series[i], tmp_trform)

    # crop images to remove black borders.
    imShape = image_registered[0].shape
    cropList = []
    for j in range(len(image_series)):
        cropped = image_registered[j][(max_cumm_shift_row+additionalCropPxl):(imShape[0]-(np.abs(min_cumm_shift_row)+additionalCropPxl)),(max_cumm_shift_col+additionalCropPxl):(imShape[1]-(np.abs(min_cumm_shift_col)+additionalCropPxl))]
        cropped = img_as_uint(cropped/np.max(cropped))
        cropList.append(cropped)
    cropList = np.stack(cropList,axis = 0)

    shiftVecDF = pd.DataFrame(shift_vectors)
    shiftVecDF.columns = ('shiftRow','shiftColumn')

    return(cropList,shiftVecDF)

def getBackground(path_to_images, im_out="./", channel="A594",sig = 5):

    files = glob.glob(path_to_images+"*"+channel+".tif")

    for imnum in range(len(files)):
        multitiff = io.MultiImage(files[imnum])
        image_series = io.collection.concatenate_images(multitiff)
        tmp = np.sum(image_series,axis=0)/np.sum(np.sum(image_series,axis=0))
        if (imnum == 0):
            outRegisters_red_sum = tmp/np.max(tmp)
        else:
            outRegisters_red = tmp/np.max(tmp)
            outRegisters_red_sum = outRegisters_red_sum + outRegisters_red

    outRegisters_red_sum = filters.gaussian(outRegisters_red_sum, sigma=(sig, sig), multichannel=False)
    outim = img_as_uint(outRegisters_red_sum/np.max(outRegisters_red_sum))

    io.imsave(im_out,outim)
    return(outim)

def correctBackground(image_series, background_image):
    background = background_image
    const_background = (background/np.sum(background))

    corrected_series = []

    for frame in range(image_series.shape[0]):
        tmpIm = np.array(np.copy(image_series[frame]))
        tmpIm = tmpIm/np.sum(tmpIm)
        tmpIm = (tmpIm-const_background)
        tmpIm[tmpIm<0] = 0
        corrected_series.append(img_as_uint(tmpIm/np.max(tmpIm)))

    corrected_series = np.stack(corrected_series,axis = 0)
    return(corrected_series)

def denoiseAndConstrast(image_series, sigMulti=15):
    patch_kw = dict(patch_size=5,
                    patch_distance=6,
                    multichannel=False)

    corrected_series = []

    for frame in range(image_series.shape[0]):
        sigma_est = np.mean(restoration.estimate_sigma(image_series[frame].astype(float), multichannel=False))
        denoise = restoration.denoise_nl_means(image_series[frame].astype(float), h=2*sigma_est, fast_mode=True,
                                **patch_kw)
        denoise = (denoise/np.max(denoise))
        denoise_contrast = exposure.equalize_adapthist(denoise)

        corrected_series.append(img_as_uint(denoise_contrast))

    corrected_series = np.stack(corrected_series,axis = 0)
    return(corrected_series)

def constrastImage(image_series):

    corrected_series = []

    for frame in range(image_series.shape[0]):
        im1 = (image_series[frame]/np.max(image_series[frame]))
        contrast = exposure.equalize_adapthist(im1)

        corrected_series.append(img_as_uint(contrast))

    corrected_series = np.stack(corrected_series,axis = 0)
    return(corrected_series)
