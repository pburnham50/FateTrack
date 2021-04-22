### Use with environment unet2DE

### script will match a corrected timelapse image to the A594-channel HCR image.



import numpy as np
from bin.fatetrack_register_v3 import *
from urllib.parse import urlparse
import cellpose
from cellpose import utils, io,models
import matplotlib
import matplotlib.pyplot as plt
import time, os, sys
import pandas as pd
import glob
### Part 1 : Image Registration
from skimage import img_as_uint,io,registration,transform,filters,restoration,util,feature,morphology,exposure,measure
from sklearn.cluster import KMeans
from scipy import ndimage
from skimage.util import montage
from scipy.spatial.distance import cdist

patch_kw = dict(patch_size=5,patch_distance=6, multichannel=False)


def getSubImage(image_path, image_name, channel='DAPI', dimx_0=0,dimx_1=1,dimy_0=0,dimy_1=1):
    hcr_best_channel = io.imread(image_path+image_name+channel+'.tif')[dimx_0:dimx_1,dimy_0:dimy_1]
    return(hcr_best_channel)

def TL_maxIntensity(nuclearImage, frame=-1):
    time_image = io.imread(nuclearImage)
    time_image = time_image[frame]
    time_image = time_image / np.max(time_image)

    tmp_time_image = morphology.white_tophat(time_image, morphology.disk(12))
    tmp_time_image = tmp_time_image/np.max(tmp_time_image)
    tmp_time_image = exposure.equalize_adapthist(tmp_time_image)
    sigma_est = np.mean(restoration.estimate_sigma(tmp_time_image.astype(float), multichannel=False))
    tmp_time_image = restoration.denoise_nl_means(tmp_time_image.astype(float), h=2*sigma_est, fast_mode=True,**patch_kw)

    time_filter = KMeans(n_clusters = 2, random_state = 0).fit_predict(tmp_time_image.reshape(np.prod(tmp_time_image.shape), 1)).reshape(tmp_time_image.shape)
    image0 = tmp_time_image*(time_filter == 0)
    image1 = tmp_time_image*(time_filter == 1)
    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    image_time_mask = time_filter == np.argmax([i0,i1])

    image_time_mask = morphology.binary_dilation(image_time_mask,morphology.diamond(1))
    image_time_mask = morphology.remove_small_objects(image_time_mask, 9)

    distance_time = ndimage.distance_transform_edt(image_time_mask)
    local_maxi_time = feature.peak_local_max(distance_time, indices=False,min_distance=9)

    io.imsave('tmp/'+nuclearImage.split('/')[-1].split('.')[0]+'_localmax.tif',local_maxi_time)

    return(local_maxi_time)

def HCR_maxIntensity(nuclearImage, saveMask=True):

    hcr_image = io.imread(nuclearImage)
    hcr_image = hcr_image / np.max(hcr_image)

    hcr_mask = 'HCRmasks/'+nuclearImage.split('/')[-1].split('.')[0]+'_locusMask.tif'

    if (os.path.exists(hcr_mask)):
        local_maxi = io.imread(hcr_mask)
    else:
        tmp_hcr_image = morphology.white_tophat(hcr_image, morphology.disk(12))
        tmp_hcr_image = tmp_hcr_image/np.max(tmp_hcr_image)
        tmp_hcr_image = exposure.equalize_adapthist(tmp_hcr_image)
        sigma_est = np.mean(restoration.estimate_sigma(tmp_hcr_image.astype(float), multichannel=False))
        tmp_hcr_image = restoration.denoise_nl_means(tmp_hcr_image.astype(float), h=2*sigma_est, fast_mode=True,**patch_kw)

        image_bkg_filter = KMeans(n_clusters = 2, random_state = 0).fit_predict(tmp_hcr_image.reshape(np.prod(tmp_hcr_image.shape), 1)).reshape(tmp_hcr_image.shape)
        image0 = tmp_hcr_image*(image_bkg_filter == 0)
        image1 = tmp_hcr_image*(image_bkg_filter == 1)
        i0 = np.average(image0[image0 > 0])
        i1 = np.average(image1[image1 > 0])
        image_bkg_filter_mask = image_bkg_filter == np.argmax([i0,i1])

        image_bkg_filter_mask = morphology.binary_dilation(image_bkg_filter_mask,morphology.diamond(1))
        image_bkg_filter_mask = morphology.remove_small_objects(image_bkg_filter_mask, 30)

        distance = ndimage.distance_transform_edt(image_bkg_filter_mask)
        local_maxi = feature.peak_local_max(distance, indices=False,min_distance=9)

        if (saveMask):
            io.imsave(hcr_mask,local_maxi)

    return(local_maxi)

def restrictSearch(TLnuclearImage,TL_maxima, HCR_maxima, totalRows, totalCols, pxOverlap):
    framename = TLnuclearImage
    frameNum = int(framename.split('/')[-1].split('_')[-3][2:])
    imDim_0,imDim_1 = np.shape(TL_maxima)
    HCRdim_0,HCRdim_1 = np.shape(HCR_maxima)

    rowNumber = int(np.ceil(frameNum/totalRows))
    colNumber = frameNum%(totalCols)
    if(colNumber==0):
        colNumber = totalCols
    if ((rowNumber%2)==0):
        colNumber = (totalCols-colNumber+1)

    pxMin_0 = int((imDim_0)*(1-pxOverlap)*(colNumber-1))
    pxMax_0 = int(HCRdim_0-((imDim_0)*(1-pxOverlap)*(totalCols-colNumber)))
    pxMin_1 = int((imDim_1)*(1-pxOverlap)*(rowNumber-1))
    pxMax_1 = int(HCRdim_1-((imDim_1)*(1-pxOverlap)*(totalRows-rowNumber)))

    return([[pxMin_0,pxMax_0],[pxMin_1,pxMax_1]])

# returns coordinates of subimage for best match of HCR to the timelapse image.
def coarseMatchHCR(TLnuclearImage,TL_maxima, HCR_maxima, stitch_tiles,penalty_multiplier,restrictedParams, cost_inclusion = 1.25):
    subdivide = stitch_tiles**2
    remove_end = stitch_tiles
    print(penalty_multiplier,stitch_tiles)
    time_positions = np.asarray(np.nonzero(TL_maxima.astype(int)))
    xcost,ycost = int(np.floor(HCR_maxima.shape[0]/subdivide)),int(np.floor(HCR_maxima.shape[1]/subdivide))
    cost = np.array(np.zeros((subdivide,subdivide)))

    newMin_0=int(np.floor((restrictedParams[0][0]/np.shape(HCR_maxima)[0])*cost.shape[0]))
    newMax_0=int(np.ceil((restrictedParams[0][1]/np.shape(HCR_maxima)[0])*cost.shape[0]))
    newMin_1=int(np.floor((restrictedParams[1][0]/np.shape(HCR_maxima)[1])*cost.shape[1]))
    newMax_1=int(np.ceil((restrictedParams[1][1]/np.shape(HCR_maxima)[1])*cost.shape[1]))

    for i in range(newMin_0,newMax_0-remove_end):
        for j in range(newMin_1,newMax_1-remove_end):
                #print(i,j)
                hcr_positions = np.asarray(np.nonzero(HCR_maxima[(j*ycost):(((j+1)*xcost)+TL_maxima.shape[1]),\
                                                                (i*xcost):(((i+1)*ycost)+TL_maxima.shape[0])]))

                costsub = np.sum(cdist(list(zip(time_positions[0], time_positions[1])),list(zip(hcr_positions[0], hcr_positions[1]))).min(axis=1))

                addedCost = ((hcr_positions.shape[1])-time_positions.shape[1])*penalty_multiplier
                cost[i,j] = costsub + addedCost

    cost[cost == 0]= np.max(cost)

    cost_list = cost.flatten()
    cost_list = np.sort(cost_list[cost_list != np.max(cost_list)])

    multiCost = (cost_list/np.min(cost_list))
    multiCost = multiCost[multiCost<=cost_inclusion]

    if (len(multiCost)>1):
        print("!!! Warning multiple points (" + str(len(multiCost)) + ") near minimum cost !!!")
    print("Cost = " + str(np.min(cost_list)/time_positions.shape[1]))

    tmpCost = cost_list[0]
    posRow = np.where(cost == tmpCost)[0][0]
    posCol = np.where(cost == tmpCost)[1][0]
    best_dims_0_0 = ((posRow*xcost)-(3*xcost))
    best_dims_0_1 = ((posRow+1)*xcost)+(3*xcost)+TL_maxima.shape[0]
    best_dims_1_0 = ((posCol*ycost)-(3*ycost))
    best_dims_1_1 = ((posCol+1)*ycost)+(3*ycost)+TL_maxima.shape[1]

    best_dims_0_0 = np.max([best_dims_0_0,0])
    best_dims_0_1 = np.min([best_dims_0_1,np.shape(HCR_maxima)[0]])
    best_dims_1_0 = np.max([best_dims_1_0,0])
    best_dims_1_1 = np.min([best_dims_1_1,np.shape(HCR_maxima)[1]])

    #hcr_positions_best = (HCR_maxima)[best_dims_0_0:best_dims_0_1,best_dims_1_0:best_dims_1_1]

    return([[best_dims_0_0,best_dims_0_1],[best_dims_1_0,best_dims_1_1]])

def fineMatchHCR(TLnuclearImage, HCRnuclearImage, coarseCoords, frame=-1):
    time_image = io.imread(TLnuclearImage)
    time_image = time_image[frame]
    time_image = time_image / np.max(time_image)

    hcr_image = io.imread(HCRnuclearImage)
    hcr_image = hcr_image / np.max(hcr_image)

    hcr_image_sub = hcr_image[coarseCoords[1][0]:coarseCoords[1][1],coarseCoords[0][0]:coarseCoords[0][1]]
    abc = np.log10(hcr_image_sub+1)

    km_filter = KMeans(n_clusters = 3, random_state = 0).fit_predict(abc.reshape(np.prod(abc.shape), 1)).reshape(abc.shape)
    image0 = abc*(km_filter == 0)
    image1 = abc*(km_filter == 1)
    image2 = abc*(km_filter == 2)

    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    i2 = np.average(image2[image2 > 0])

    hcr_image_sub = hcr_image_sub*(km_filter!=np.argmin([i0,i1,i2]))

    pad_x = (hcr_image_sub.shape[0] - time_image.shape[0])
    pad_y = (hcr_image_sub.shape[1] - time_image.shape[1])
    time_image_pad = np.pad(time_image,((pad_x,0),(pad_y,0)),mode='constant')

    crossCor = registration.phase_cross_correlation(time_image_pad,hcr_image_sub)
    #tmp_trform = transform.AffineTransform(translation=[-crossCor[0][0],-crossCor[0][1]])
    #image_registered = transform.warp(hcr_image_sub, tmp_trform)
    #image_registered = image_registered[pad_x:,pad_y:]

    translate_dims_0_0 = (coarseCoords[0][0]-int(crossCor[0][1]))+pad_x
    translate_dims_0_1 = (coarseCoords[0][1]-int(crossCor[0][1]))
    translate_dims_1_0 = (coarseCoords[1][0]-int(crossCor[0][0]))+pad_y
    translate_dims_1_1 = (coarseCoords[1][1]-int(crossCor[0][0]))

    #hcr_fine_match = hcr_image[translate_dims_0_0:translate_dims_0_1,translate_dims_1_0:translate_dims_1_1]

    #hcr_fine_match = exposure.equalize_adapthist(hcr_fine_match)

    return([[translate_dims_0_0,translate_dims_0_1],[translate_dims_1_0,translate_dims_1_1]])

def getHCRsubImage(outFile, TLnuclearImage, HCRnuclearImage, Fine_HCRsubcoordsFile, Coarse_HCRsubcoordsFile, Restrict_HCRsubcoordsFile,stitch_tiles,penalty_multiplier,cost_inclusion, totalRows, totalCols,pixelOverlap,frame=-1):
    hcr_loci = HCR_maxIntensity(nuclearImage=HCRnuclearImage, saveMask=True)
    tl_loci = TL_maxIntensity(nuclearImage=TLnuclearImage, frame=frame)
    restrictCoords = restrictSearch(TLnuclearImage=TLnuclearImage,TL_maxima=tl_loci, HCR_maxima=hcr_loci, totalRows=totalRows, totalCols=totalCols, pxOverlap=pixelOverlap)
    coarseCoords = coarseMatchHCR(TLnuclearImage=TLnuclearImage,TL_maxima=tl_loci, HCR_maxima=hcr_loci, stitch_tiles=stitch_tiles,penalty_multiplier=penalty_multiplier,cost_inclusion = cost_inclusion, restrictedParams=restrictCoords)
    fineMatchcoords = fineMatchHCR(TLnuclearImage=TLnuclearImage, HCRnuclearImage=HCRnuclearImage, coarseCoords=coarseCoords, frame=frame)
    hcr_image = io.imread(HCRnuclearImage)
    hcr_image = hcr_image / np.max(hcr_image)
    fineMatch = hcr_image[fineMatchcoords[1][0]:fineMatchcoords[1][1],fineMatchcoords[0][0]:fineMatchcoords[0][1]]
    fineMatch = exposure.equalize_adapthist(fineMatch)
    #fineMatch = fineMatch.astype('uint8')
    fineMatch = img_as_uint(fineMatch / np.max(fineMatch))
    io.imsave(outFile,fineMatch)

    np.savetxt(Restrict_HCRsubcoordsFile, np.array(restrictCoords), fmt="%s")
    np.savetxt(Coarse_HCRsubcoordsFile, np.array(coarseCoords), fmt="%s")
    np.savetxt(Fine_HCRsubcoordsFile, np.array(fineMatchcoords), fmt="%s")
    return(coarseCoords)

def getHCRmask(HCRsubImage, nucDiameter=0., cellprob_threshold=0., flow_threshold=.4):

    model = models.Cellpose(gpu=False, model_type='nuclei')
    image = io.imread(HCRsubImage)
    image = exposure.equalize_adapthist(image)

    masks, flows, styles, diams = model.eval(image, diameter=nucDiameter, \
                                            cellprob_threshold=cellprob_threshold, \
                                            flow_threshold=flow_threshold, channels=[0,0])

    cellpose.io.masks_flows_to_seg(images=image, masks=masks, flows=flows, diams=diams, channels=[0,0], file_names=HCRsubImage)
