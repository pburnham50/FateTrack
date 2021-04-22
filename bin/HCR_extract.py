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
from scipy import ndimage as ndi
from skimage.util import montage
from skimage.transform import warp_polar, rotate, rescale
from scipy.spatial import Voronoi
from skimage.feature import peak_local_max

from sklearn.cluster import KMeans

matplotlib.use('agg')

def getSubImage(hcr_file, hcr_coords_file, channel='DAPI'):
    image_name='_'.join('/'.join(hcr_file.split('.')[:-1]).split('_')[:-1])+'_'
    HCRcoords = np.loadtxt(hcr_coords_file).astype(int)
    hcr_best_channel = io.imread(image_name+channel+'.tif')[HCRcoords[0]:HCRcoords[1],HCRcoords[2]:HCRcoords[3]].astype('float32')
    hcr_best_channel = rotate(hcr_best_channel, HCRcoords[4], mode='constant').astype('uint16')
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
    temp = np.asarray(np.load(HCR_mask_file,allow_pickle=True))#.item()
    masks = temp.T #segmentation.clear_border(temp)#['masks'])
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


def getVoronoiStyle(seg_file,max_voro_area,voro_imfile,voro_imfile_2,voro_outfile,voro_transfile):
    hcr_file = '/'.join(seg_file.split('/')[:-1])+'/HCRsub_A594.tif'
    im = io.imread(hcr_file)
    im = np.zeros_like(np.array(im))

    temp = np.asarray(np.load(seg_file,allow_pickle=True)).item()
    masks = temp['masks']

    fro = pd.DataFrame(measure.regionprops_table(masks, properties=['label','centroid']))


    points_mask = np.array(fro[['centroid-0','centroid-1']].to_numpy())

    vor = Voronoi(points_mask)

    my_dpi=im.shape[1]

    plt.rcParams['figure.dpi'] = my_dpi
    plt.rcParams['figure.figsize'] = ( im.shape[0]/my_dpi,im.shape[1]/my_dpi)
    fig = plt.figure();

    for simplex in vor.ridge_vertices:
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-',c='black',linewidth=.2)

    center = points_mask.mean(axis=0)
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.any(simplex < 0):
            i = simplex[simplex >= 0][0] # finite end Voronoi vertex
            t = points_mask[pointidx[0]] - points_mask[pointidx[1]]  # tangent
            t = t / np.linalg.norm(t)
            n = np.array([-t[1], t[0]]) # normal
            midpoint = points_mask[pointidx].mean(axis=0)
            far_point = vor.vertices[i] + np.sign(np.dot(midpoint - center, n)) * n * 100
            plt.plot([vor.vertices[i,0], far_point[0]],
                     [vor.vertices[i,1], far_point[1]], 'k-',c='black',linewidth=.2)

    plt.xlim([0, im.shape[0]]); plt.ylim([0,im.shape[1]])
    plt.axis('off')
    fig.tight_layout(pad=0)
    plt.savefig(voro_imfile, dpi=my_dpi, #bbox_inches='tight',#dpi=my_dpi,
                transparent=False, pad_inches=0,facecolor='white')
    plt.close()
    im2 = io.imread(voro_imfile)
    voro = (im2[:,:,0])
    voro = voro[1:-1, 1:-1]
    voro = np.pad(voro, pad_width=1, mode='constant')

    distance = ndi.distance_transform_edt(voro)
    coords = peak_local_max(distance, footprint=np.ones((1, 1)), labels=voro)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(mask)
    labels = segmentation.watershed(-distance, markers, mask=voro)
    labels = morphology.remove_small_objects(labels, min_size=40, connectivity=1, in_place=False)
    labels = morphology.dilation(labels, morphology.square(3))

    segmasks = masks
    segmasks = morphology.dilation(segmasks,morphology.square(3))

    sizeOfSegs = pd.DataFrame(measure.regionprops_table(labels, properties=['label','area']))
    bigMasks = np.array(sizeOfSegs[sizeOfSegs['area']>=max_voro_area]['label'])

    newVorMask = np.copy(labels)[::-1,:]
    for bMI in range(len(bigMasks)):

        chckMtx = (labels == bigMasks[bMI])[::-1,:]

        for i in range(len(points_mask)):
            if(chckMtx[int(np.round(points_mask[i][1])),int(np.round(points_mask[i][0]))]):
                confirm = points_mask[i]

        tmp_cellpose_mask = (morphology.dilation((segmasks == int(fro[(fro['centroid-0']==confirm[0])&(fro['centroid-1']==confirm[1])]['label'])).T,morphology.disk(11))).astype(int)
        tmp_voronoi_mask = 2*chckMtx.astype(int)
        tmp_join = segmentation.join_segmentations(tmp_cellpose_mask,tmp_voronoi_mask)
        tmp_join = (tmp_join == np.max(tmp_join))

        newVorMask[newVorMask == bigMasks[bMI]] = 0
        newVorMask[tmp_join] = bigMasks[bMI]

    np.save(voro_outfile, newVorMask, allow_pickle=True, fix_imports=True)
    io.imsave(voro_imfile_2, segmentation.find_boundaries(newVorMask))

    oldAssign = pd.DataFrame(measure.regionprops_table(masks, properties=['label','centroid']))
    newAssign = pd.DataFrame(measure.regionprops_table(newVorMask, properties=['label','centroid']))

    Clps2Voro = pd.DataFrame()

    for nlab in range(newAssign.shape[0]):
        tmpMtx = (newVorMask == newAssign['label'][nlab])
        for olab in range(oldAssign.shape[0]):
            if (tmpMtx[int(np.round(oldAssign['centroid-1'][olab])),int(np.round(oldAssign['centroid-0'][olab]))]):
                Clps2Voro = Clps2Voro.append(pd.DataFrame([newAssign['label'][nlab], oldAssign['label'][olab]]).T)

    Clps2Voro = Clps2Voro.rename(columns={0: "voro_label", 1: "clps_label"})
    Clps2Voro = Clps2Voro.reset_index(drop=True)
    Clps2Voro.to_csv(voro_transfile)


def assembleFinalState(HCR_mask_file, HCR_image_file, HCR_coords_file, Voro_mask_file,Voro_image_file,Voro_image_file_final, Voro_Transfer_file,channelList, final_timepoint,dilations=[1],MaxVoroArea = 1600):

    result_subimages = list(map(lambda x : getSubImage(HCR_image_file, hcr_coords_file=HCR_coords_file, channel=x), channelList))
    result_bkgrd = list(map(getBackgroundValue, result_subimages))

    fullExpChannelFrame = pd.DataFrame()

    getVoronoiStyle(seg_file=HCR_mask_file,max_voro_area=MaxVoroArea,
                    voro_imfile=Voro_image_file,voro_imfile_2 =Voro_image_file_final,
                    voro_outfile=Voro_mask_file,voro_transfile=Voro_Transfer_file)

    for chan in range(len(channelList)):
        for dil in dilations:
            if (len(fullExpChannelFrame)==0):
                fullExpChannelFrame = getIntProps(HCR_mask_file = Voro_mask_file,
                                                  channel = channelList[chan],
                                                  channelImageArray = result_subimages[chan],
                                                  dilation=dil, bkgrdValue=result_bkgrd[chan])
            else:
                tmp_fullExpChannelFrame = getIntProps(HCR_mask_file = Voro_mask_file,
                                                  channel = channelList[chan],
                                                  channelImageArray = result_subimages[chan],
                                                  dilation=dil, bkgrdValue=result_bkgrd[chan])
                fullExpChannelFrame = pd.merge(fullExpChannelFrame,tmp_fullExpChannelFrame)
    fullExpChannelFrame['Voro_ID_'+str(final_timepoint)] = str(final_timepoint)+'_'+fullExpChannelFrame['label'].astype(str)
    return(fullExpChannelFrame)

def getAllConnections(pathToConnects, finalConnection, final_timepoint, sampleName):
    connections = pd.read_csv(pathToConnects+'/'+'Connections_MasterID_'+str(0)+'.csv')
    for imnum in range(1,final_timepoint-1):
        connections = pd.merge(connections,pd.read_csv(pathToConnects+'/'+'Connections_MasterID_'+str(imnum)+'.csv'))
    connections['TrackID'] = sampleName+':'+connections.agg('.'.join, axis=1)
    finalConnect = pd.read_csv(finalConnection)
    connections = pd.merge(connections,finalConnect)

    return(connections)
