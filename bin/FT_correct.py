"""
Philip Burnham
Arjun Raj Lab
FateTrack - Image Correction and Registration
2021
FT_correct.py
"""

import argparse
import matplotlib
from fatetrack_register_v3 import *
import time, os, sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
from skimage import io,registration,transform,filters,img_as_uint,restoration


def main():
    parser = argparse.ArgumentParser('Register images using phase cross correlation.')

    parser.add_argument('--nucImage', type = str, default='./',dest='nucImage',help = 'Path and name of nuclear-fluorescent image.')
    parser.add_argument('--bfImage', type = str, default='./',dest='bfImage',help = 'Path and name of brightfield image.')
    parser.add_argument('--outpath', type = str, default='./',dest='outpath',help = 'Path to write out images.')
    parser.add_argument('--backgroundImage', type = str, dest='background', help = 'None if no correction, put directory of group of images to calculate, put image path and filename for single image.')
    parser.add_argument('--cropPlus', type = int, default=2, dest='cropPlus',help = 'Additional pixel border to crop.')
    parser.add_argument('--startFrame', type = int, default=1, dest='startFrame',help = 'Start frame number.')
    parser.add_argument('--stopFrame', type = int, default=1000, dest='stopFrame',help = 'Final frame number.')
    args = parser.parse_args()

    backgroundImage = args.background
    nucImage = args.nucImage
    bfImage = args.bfImage
    outpath = args.outpath
    pxCrop = args.cropPlus
    startFrame = args.startFrame
    stopFrame = args.stopFrame

    fov = '_'.join(nucImage.split('/')[-1].split('.')[0].split('_')[:-1])

    multitiff = io.MultiImage(nucImage)
    image_series = io.collection.concatenate_images(multitiff)

    multitiff_bf = io.MultiImage(bfImage)
    image_series_bf = io.collection.concatenate_images(multitiff_bf)

    if ((startFrame>1)|(stopFrame<(image_series.shape[0]))):
        image_series = image_series[range(startFrame-1,stopFrame)]
        image_series_bf = image_series_bf[range(startFrame-1,stopFrame)]

    if (backgroundImage == "None"):
        corrected_series = image_series
    elif (backgroundImage[-1] == "/"):
        bkgrd = getBackground(backgroundImage, "./"+fov+"_A594_background.tif",channel="A594",sig = 1)
        corrected_series = correctBackground(image_series, bkgrd)
    else:
        bkgrd = io.imread(backgroundImage)
        corrected_series = correctBackground(image_series, bkgrd)


    bf_registered_cells = register(image_series=image_series_bf, additionalCropPxl=pxCrop)
    red_registered_cells = applyRegistration(image_series=corrected_series, shiftFrame=bf_registered_cells[1], additionalCropPxl=pxCrop)

    denoise_red = denoiseAndConstrast(red_registered_cells[0], sigMulti=15)
    denoise_bf = denoiseAndConstrast(bf_registered_cells[0], sigMulti=15)

    if  1-(os.path.isdir(outpath)):
            os.mkdir(outpath)

    io.imsave(outpath + fov + "_Brightfield_corr.tif",denoise_bf)
    io.imsave(outpath + fov + "_Nuclear_corr.tif",denoise_red)

    return

if __name__ == '__main__':
    main()
