### HCR_subImage.smk
from bin.HCR_matching import *
from skimage import io
import cv2
import numpy as np
from skimage.registration import phase_cross_correlation
from skimage.transform import warp_polar, rotate, rescale

# Connect frames

rule HCR_subImage:
  input:
    TLimage='{path}/correction/{sample}_Nuclear_corr.tif',
    HCRimage=config['experiment']['HCR_A594']
  output:
    HCRsub='{path}/HCR/{sample}_finalState_Nuclear.tif',
    HCRcoords='{path}/HCR/{sample}_HCRsubCoords.txt',
    HCRmask='{path}/HCR/{sample}_finalState_Nuclear_seg.npy'
  params:
    endFrame=config['correction']['endFrame'],
    flowThreshold=config['segment']['flowThreshold'],
    cellProbabilityThreshold=config['segment']['cellProbabilityThreshold'],

    nucDiam=config['segment']['nuclearDiameter']
  run:
    hcr_image = io.imread(input.HCRimage)
    lastTL = io.imread(input.TLimage)[(params.endFrame - 1)]
    hcr_image = hcr_image.astype('float32')
    lastTL = lastTL.astype('float32')

    result = cv2.matchTemplate(hcr_image,lastTL,cv2.TM_CCOEFF_NORMED)
    point = np.unravel_index(result.argmax(),result.shape)

    coord00,coord01 = point[0],(point[0]+lastTL.shape[0])
    coord10,coord11 = point[1],(point[1]+lastTL.shape[1])

    HCRsubimage = hcr_image[coord00:coord01,coord10:coord11]

    TL_polar = warp_polar(lastTL)
    HCR_polar = warp_polar(HCRsubimage)

    shifts, error, phasediff = phase_cross_correlation(HCR_polar, TL_polar)

    np.savetxt(output.HCRcoords, np.array([coord00,coord01,coord10,coord11, shifts[0]]), fmt="%s")
    io.imsave(output.HCRsub,rotate(HCRsubimage, shifts[0]).astype('uint16'))

    getHCRmask(HCRsubImage=output.HCRsub, nucDiameter=params.nucDiam, \
                cellprob_threshold=params.cellProbabilityThreshold, \
                flow_threshold=params.flowThreshold) ;
