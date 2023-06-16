"""
Philip Burnham
Arjun Raj Lab
FateTrack - nuclear feature extraction
2021
fatetrack_nucFeatureExtraction.py
"""
# Import various libraries
import numpy as np
import time, os, sys, math
import matplotlib.pyplot as plt
import glob
from scipy import ndimage as ndi
from skimage import color, feature, filters, io, measure, morphology, restoration, segmentation, exposure, restoration, data, exposure
import warnings
import pandas as pd
from skimage.feature import greycomatrix, greycoprops
from numpy import genfromtxt
from skimage.segmentation import clear_border
import mahotas
from skimage.feature import hog

featureColumns = ('label','area','filled_area','bbox','bbox_area', 'centroid','convex_area','eccentricity',
                      'solidity','convex_area','mean_intensity','min_intensity','filled_area',
                      'max_intensity','orientation','major_axis_length','minor_axis_length','euler_number',
                      'perimeter','extent','intensity_image','equivalent_diameter',
                      'moments_central','moments_hu','moments','moments_normalized',
                      'inertia_tensor','inertia_tensor_eigvals')

def extraObjectTraits(intensityImage):

    im_stretched = mahotas.stretch(intensityImage)

    zernike = mahotas.features.zernike(im_stretched, radius=20, degree=8, cm=mahotas.center_of_mass(im_stretched))
    zern_df = pd.DataFrame(zernike).T
    zern_df.columns = ['zern_'+ sub for sub in list(map(str,zern_df.columns))]

    haralick = mahotas.features.haralick(im_stretched,return_mean=True)
    hara_df = pd.DataFrame(haralick).T
    hara_df.columns = ['hara_'+ sub for sub in list(map(str,hara_df.columns))]
    width,height = np.shape(im_stretched)
    maxsize=80
    h1 = int(np.floor((maxsize-height)/2))
    h2 = int(np.ceil((maxsize-height)/2))
    w1 = int(np.floor((maxsize-width)/2))
    w2 = int(np.ceil((maxsize-width)/2))
    padIm =(np.pad(im_stretched,pad_width=((w1,w2),(h1,h2)),mode='constant'))

    hog_fd = hog(padIm, orientations=20, pixels_per_cell=(25, 25), cells_per_block=(2, 2), visualize=False, multichannel=False,feature_vector=True)
    hog_df = pd.DataFrame(hog_fd).T
    hog_df.columns = ['hog_'+ sub for sub in list(map(str,hog_df.columns))]

    return(pd.concat([hog_df, hara_df,zern_df], axis=1))

def collectObjectTraits(segmentationsMask, channelImage):

    image_props = pd.DataFrame(measure.regionprops_table(segmentationsMask, intensity_image=channelImage,properties=featureColumns))
    im_df = pd.DataFrame(image_props)

    additionProps = pd.DataFrame()

    for object in range(im_df.shape[0]):
        tmpProps = extraObjectTraits(im_df['intensity_image'][object])
        additionProps = additionProps.append(tmpProps)

    additionProps = additionProps.reset_index(drop=True)
    final = pd.concat([im_df, additionProps], axis=1)
    return(final)

def staticFeatureExt(pathToSegmenations, imageFile, outpath, sampleName):

    if  1-(os.path.isdir(outpath)):
                os.mkdir(outpath)

    numFrames = len(glob.glob(pathToSegmenations+'/*_seg.npy'))
    imorig = io.imread(imageFile);

    for frame in range(numFrames):
        segmentArrays = np.asarray(np.load(pathToSegmenations+'/'+str(frame)+'_seg.npy',allow_pickle=True)).item()
        tmpMask = clear_border(segmentArrays['masks'])
        tmpImage = imorig[frame]
        featureFrame = collectObjectTraits(segmentationsMask=tmpMask, channelImage=tmpImage)
        featureFrame["frame"] = frame

        export = featureFrame.to_csv(outpath + sampleName +'.'+ str(frame) + "_staticFeatures.csv", index = None, header=True)
