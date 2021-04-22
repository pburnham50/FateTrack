### testing image registration functions
import numpy as np
from urllib.parse import urlparse
from cellpose import utils, io,models
import matplotlib
import matplotlib.pyplot as plt
import time, os, sys
import pandas as pd
import glob
### Part 1 : Image Registration
from skimage import img_as_uint,io,registration,transform,filters,restoration,util,feature,morphology,exposure,measure,segmentation
from sklearn.cluster import KMeans
from scipy import ndimage
from skimage.util import montage

from sklearn.cluster import KMeans

def getSubImage(hcr_file, hcr_coords_file, channel='DAPI'):
    image_name='_'.join('/'.join(hcr_file.split('.')[:-1]).split('_')[:-1])+'_'
    HCRcoords = np.loadtxt(hcr_coords_file).astype(int)
    hcr_best_channel = io.imread(image_name+channel+'.tif')[HCRcoords[1][0]:HCRcoords[1][1],HCRcoords[0][0]:HCRcoords[0][1]]
    io.imsave('/'.join(hcr_coords_file.split('/')[:-1])+'/HCRsub_'+channel +'.tif',hcr_best_channel)
    return(hcr_best_channel)

def getBackgroundValue(subImage):
    #abc = np.log10(subImage+1)
    #abc= np.nan_to_num(abc)
    abc=subImage
    km_filter = KMeans(n_clusters = 2, random_state = 0).fit_predict(abc.reshape(np.prod(abc.shape), 1)).reshape(abc.shape)
    image0 = abc*(km_filter == 0)
    image1 = abc*(km_filter == 1)

    i0 = np.average(image0[image0 > 0])
    i1 = np.average(image1[image1 > 0])
    return((np.min([i0,i1])).astype('int'))

def getIntProps(HCR_mask_file, channel, channelImageArray, dilation=1, bkgrdValue=0):
    temp = np.asarray(np.load(HCR_mask_file,allow_pickle=True)).item()
    masks = segmentation.clear_border(temp['masks'])
    masks = morphology.dilation(masks,morphology.square(dilation))
    props = measure.regionprops_table(masks, intensity_image=channelImageArray, properties=('label','centroid',
                                                                                            'filled_area',
                                                     'min_intensity',
                                                     'mean_intensity',
                                                     'max_intensity'))
    props = pd.DataFrame(props)
    if(dilation>1):
        props = props.drop(['centroid-0', 'centroid-1'], axis=1)
    props['sum_intensity'] =  props['filled_area'] * props['mean_intensity']
    props = props.add_suffix('_'+str(dilation))
    props = props.add_prefix(channel+'_')
    props[channel+'_background'] = bkgrdValue
    props = props.rename(columns={channel+"_label_"+str(dilation): "label",
                                  channel+"_centroid-0_"+str(dilation): "centroid-0",
                                  channel+"_centroid-1_"+str(dilation): "centroid-1",
                                 })
    return(props)


def assembleFinalState(HCR_mask_file, HCR_image_file, HCR_coords_file, channelList, final_timepoint,dilations=[1]):

    result_subimages = list(map(lambda x : getSubImage(HCR_image_file, hcr_coords_file=HCR_coords_file, channel=x), channelList))
    result_bkgrd = list(map(getBackgroundValue, result_subimages))

    fullExpChannelFrame = pd.DataFrame()

    for chan in range(len(channelList)):
        for dil in dilations:
            if (len(fullExpChannelFrame)==0):
                fullExpChannelFrame = getIntProps(HCR_mask_file = HCR_mask_file,
                                                  channel = channelList[chan],
                                                  channelImageArray = result_subimages[chan],
                                                  dilation=dil, bkgrdValue=result_bkgrd[chan])
            else:
                tmp_fullExpChannelFrame = getIntProps(HCR_mask_file = HCR_mask_file,
                                                  channel = channelList[chan],
                                                  channelImageArray = result_subimages[chan],
                                                  dilation=dil, bkgrdValue=result_bkgrd[chan])
                fullExpChannelFrame = pd.merge(fullExpChannelFrame,tmp_fullExpChannelFrame)
    fullExpChannelFrame['Master_ID_'+str(final_timepoint)] = str(final_timepoint)+'_'+fullExpChannelFrame['label'].astype(str)
    return(fullExpChannelFrame)

def getAllConnections(pathToConnects, finalConnection, final_timepoint, sampleName):
    connections = pd.read_csv(pathToConnects+'/'+'Connections_MasterID_'+str(0)+'.csv')
    for imnum in range(1,final_timepoint-1):
        connections = pd.merge(connections,pd.read_csv(pathToConnects+'/'+'Connections_MasterID_'+str(imnum)+'.csv'))
    connections['TrackID'] = sampleName+':'+connections.agg('.'.join, axis=1)
    finalConnect = pd.read_csv(finalConnection)
    connections = pd.merge(connections,finalConnect)

    return(connections)
