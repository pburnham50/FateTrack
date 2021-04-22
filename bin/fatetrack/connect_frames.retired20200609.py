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


def buildFeatureFrame(filename,timepoint):
    temp = np.asarray(np.load(filename,allow_pickle=True)).item()
    image_props = measure.regionprops_table(temp['masks'],
                                        intensity_image=temp['img'],
                                        properties=('label','area','filled_area','bbox', 'centroid',
                                                    'eccentricity','solidity','convex_area',
                                                    'mean_intensity','min_intensity','max_intensity',
                                                    'orientation','major_axis_length','minor_axis_length',
                                                    'perimeter','extent','intensity_image'))
    im_df = pd.DataFrame(image_props)
    im_df['time'] = timepoint
    return(im_df)

def buildOffsetFrame(filename_tminus1,filename_tplus1):
    temp1 = np.asarray(np.load(filename_tminus1,allow_pickle=True)).item()
    temp2 = np.asarray(np.load(filename_tplus1,allow_pickle=True)).item()

    image_props = measure.regionprops_table(temp1['masks'],
                                        intensity_image=temp2['img'],
                                        properties=('label','centroid','area',"mean_intensity"))
    im_df = pd.DataFrame(image_props)
    im_df['time'] = None
    return(im_df)

def linkEnergy(image1, image2, im1_select, im2_select):
    deltaX = np.sqrt((image1['centroid-0'][im1_select]-image2['centroid-0'][im2_select])**2+
                     (image1['centroid-1'][im1_select]-image2['centroid-1'][im2_select])**2)
    deltaA = np.absolute(image1['area'][im1_select] - image2['area'][im2_select])
    score = deltaX + np.sqrt(deltaA)
    return(score)

def generateLinks(image1, image2, im1_select, dist_multiplier=2):
    delX = np.sqrt((image1['centroid-0'][im1_select]-image2['centroid-0'])**2+
                   (image1['centroid-1'][im1_select]-image2['centroid-1'])**2)
    max_dist = dist_multiplier*min(delX)
    candidates = np.array(delX[delX < max_dist].index)
    return(candidates)

def ScoreMotherArea(motherArea, ProgArea_1, Prog_Area_2, threshold=0.95):
    if (motherArea/(ProgArea_1 + Prog_Area_2)) > threshold :
        return(0)
    else:
        return(10)

def ScoreMotherInt(motherCurrentInt, motherNextInt):
    normInt = (motherCurrentInt/(motherNextInt + 10**-4))
    if (normInt > 2):
        return(-1)
    elif (normInt <= 2)&(normInt > 1.4):
        return(0)
    else:
        return(1)

def ScoreProjInt(projCurrentInt, projPrevInt):
    normInt = (projCurrentInt/(projPrevInt + 10**-4))
    if (normInt > 2):
        return(-1)
    elif (normInt <= 2)&(normInt > 1):
        return(0)
    else:
        return(1)

def ScoreProjDiffInt(proj1Int, proj2Int):
    return(np.abs(proj1Int - proj2Int)/2)

def ScoreDivisionTotal(motherFrameCurr, motherFrameNext,
                       projFrameCurr, projFramePrev,
                        motherCell,projCell_1, projCell_2):

    motherArea = ScoreMotherArea(motherFrameCurr["area"][motherCell],
                                 projFrameCurr["area"][projCell_1],
                                 projFrameCurr["area"][projCell_2])

    motherInt = ScoreMotherInt(motherFrameCurr["mean_intensity"][motherCell],
                               motherFrameNext["mean_intensity"][motherCell])

    projInt = -1 + ScoreProjInt(projFrameCurr["mean_intensity"][projCell_1],  projFramePrev["mean_intensity"][projCell_1]) +ScoreProjInt(projFrameCurr["mean_intensity"][projCell_2], projFramePrev["mean_intensity"][projCell_2])

    projIntDiff = ScoreProjDiffInt(projFrameCurr["mean_intensity"][projCell_1],
                                   projFrameCurr["mean_intensity"][projCell_2])

    return(motherArea + motherInt + projInt + projIntDiff )

def DivisionCandidate(motherFrameCurr, motherFrameNext,
                       projFrameCurr, projFramePrev,
                        motherCell,projCell_1, projCell_2_candidates, threshold=3):
    tru_vec=[]
    for i in projCell_2_candidates:
        if(ScoreDivisionTotal(motherFrameCurr,motherFrameNext,
                              projFrameCurr,projFramePrev,
                              motherCell,projCell_1,i) < threshold):
            tru_vec=np.append(tru_vec,True)
        else:
            tru_vec=np.append(tru_vec,False)
    return(np.any(tru_vec))

def buildConnections(filename_t0,filename_t1,time0,time1,greedy=True,openingCost=2, nnDist=3, DivScoreThreshold=3):
    ip0 = buildFeatureFrame(filename_t0,time0)
    ip1 = buildFeatureFrame(filename_t1,time1)
    fx0 = buildOffsetFrame(filename_t0,filename_t1)
    fx1 = buildOffsetFrame(filename_t1,filename_t0)


    num=0
    arr = pd.DataFrame([]).T
    for i in np.array(ip0.index):
        candidates = generateLinks(ip0, ip1, i, nnDist)
        for j in range(len(candidates)):
            proj1 = candidates[j]
            proj2pairs = np.delete(candidates,j)
            if(len(proj2pairs)>0):
                divscore = (DivisionCandidate(ip0,fx0,ip1,fx1,i,candidates[j],proj2pairs,DivScoreThreshold))
            else:
                divscore = False
            arr = arr.append(pd.DataFrame([num,i,proj1,linkEnergy(ip0, ip1, i, proj1),divscore]).T)
            num=num+1

    arr.columns = ['index','prev','next','score','divisionCandidate']
    arr.index = arr['index']
    arr.iloc[np.array(arr[arr['divisionCandidate']].index),3] = np.array(arr[arr['divisionCandidate']]['score']/2)


    if(greedy ==True):
        nextFeatList = np.unique(arr['next'])
        nextCrop = pd.DataFrame()

        for next in nextFeatList:
            subarr = arr[arr['next'] == next]
            new = subarr[subarr.score == subarr.score.min()]
            nextCrop = nextCrop.append(new)

        prevFeatList = np.unique(nextCrop['prev'])
        prevCrop = pd.DataFrame()

        for prev in prevFeatList:
            subarr = nextCrop[nextCrop['prev'] == prev]
            if np.sum((subarr['divisionCandidate']).astype(int))>1 :
                new=subarr
            else:
                new = subarr[subarr.score == subarr.score.min()]
            prevCrop =  prevCrop.append(new)

        final = prevCrop
    else:
        final = arr
    final['timePrev']=time0
    final['timeNext']=time1
    final['prev'] = final['prev']+1
    final['next'] = final['next']+1
    final['PrevID'] = final['timePrev'].astype(str)+'-'+final['prev'].astype(str)
    final['NextID'] = final['timeNext'].astype(str)+'-'+final['next'].astype(str)
    final['ConnectID'] = final['PrevID']+':'+final['NextID']

    final = final.drop(['index','timePrev', 'timeNext','PrevID',"NextID"], axis=1)
    final['divisionCandidate'] = final['divisionCandidate'].astype(int)
    final = final.rename(columns={"prev": time0, "next": time1,
                                  "score": "cost", "divisionCandidate": "divisions"})

    for i in (np.unique(ip0['label'])):
        if len(final[final[str(time0)]==i])<1 :
                addin = pd.DataFrame([i,0,openingCost,0,
                                      str(time0)+"-"+str(i)+":"+str(time0)+"-0"]).T
                addin.columns = final.columns
                final = final.append(addin)

    for i in (np.unique(ip1['label'])):
        if len(final[final[str(time1)]==i])<1 :
                addin = pd.DataFrame([0,i,openingCost,0,
                                      str(time0)+"-0:"+str(time1)+"-"+str(i)]).T
                addin.columns = final.columns
                final = final.append(addin)

    addin = pd.DataFrame([0,0,1.25*openingCost,0, str(time0)+"-0:"+str(time1)+"-0"]).T
    addin.columns = final.columns
    final = final.append(addin)
    final['lastDiv'] = 0

    final = final[final.cost < (1.25*openingCost)+1]

    return(final)

def mergeConnections(ConnectFrame1, ConnectFrame2, timePoint):
    result = pd.merge(ConnectFrame1, ConnectFrame2, on=timePoint)
    result['cost'] = result['cost_x'] + result['cost_y']
    result['divisions'] = result['divisions_x'] + result['divisions_y']
    result['lastDiv'] = (1-result['divisions_y'])*result['lastDiv_x'] + result['divisions_y']*int(timePoint)
    nextPoint = str(int(timePoint)+1)
    result['ConnectID'] = result['ConnectID_x']+':' + nextPoint+"-"+result[nextPoint].astype(str)
    result = result.drop(['cost_x','cost_y', 'divisions_x','divisions_y',
                          'ConnectID_x','ConnectID_y','lastDiv_x','lastDiv_y'], axis=1)
    return(result)

def callFeatureMtx(dataFrame, frameNum, featureNum):
    return(dataFrame[frameNum][dataFrame[frameNum]['label']==featureNum])

def distanceFrom(Frame,origTime=0):
    valvec = np.sqrt((Frame['centroid-0']-Frame['centroid-0'][origTime])**2+(Frame['centroid-1']-Frame['centroid-1'][origTime])**2)
    return(valvec)

def distanceFromPrevious(Frame,TimeBefore=1):
    dpvec=[0.]
    for i in range(1,len(Frame)):
        valvec = [np.sqrt((Frame['centroid-0'][i]-Frame['centroid-0'][i-TimeBefore])**2+(Frame['centroid-1'][i]-Frame['centroid-1'][i-TimeBefore])**2)]
        dpvec = np.concatenate((dpvec,valvec))
    return(dpvec)

def displacementFrom(distanceVector):
    return(np.cumsum(distanceVector))

def callFeatureHistory(ConnectionID):
    frfets = ConnectionID.split(':')
    featDF = pd.DataFrame()
    for i in range(len(frfets)):
        featDF=featDF.append(callFeatureMtx(Frames,
                                            int(frfets[i].split('-')[0]),
                                            int(frfets[i].split('-')[1])))
    featDF.index = featDF.time

    featDF["totalDisplacement"] = distanceFrom(featDF,origTime=0)
    featDF["previousDistance"] = distanceFromPrevious(featDF,1)
    featDF["totalDistance"] = displacementFrom(featDF["previousDistance"])

    featDF["wandering"] = featDF["totalDistance"]/featDF["totalDisplacement"]
    featDF["ConnectID"] = ConnectionID

    return(featDF)

def selectMinCost(originalFrame, frameNumber):
    revised = pd.DataFrame()
    for i in (np.unique(originalFrame[frameNumber]))[1:]:
        subdf = originalFrame[originalFrame[frameNumber] == i]
        subdf = subdf[(subdf.iloc[:,1:int(frameNumber)]>0).all(axis='columns')]
        subdf = subdf[subdf['0'] != 0]

        minCost = subdf[subdf.cost == subdf.cost.min()]
        revised = revised.append(minCost)
    return(revised)


###############################################################################################################
# main function
###############################################################################################################

def main():
    parser = argparse.ArgumentParser('Generate file of image properties')

    parser.add_argument('--i1', type = str, dest='i1', help = 'Input results of cellpose for image 1.')
    parser.add_argument('--i2', type = str, dest='i2', help = 'Input results of cellpose for image 2.')
    parser.add_argument('--t1', type = str, default='0', dest='timepoint1',help = ' timepoint labelling image 1.')
    parser.add_argument('--t2', type = str, default='1', dest='timepoint2',help = ' timepoint labelling image 2.')
    parser.add_argument('--greedy', type = bool, default=True, dest='greed',help = 'Run greedy algorithm')
    parser.add_argument('--openingCost', type = int, default=25,dest='openCost',help = 'Cost for opening/closing a track.')
    parser.add_argument('--maxCost', type = int, default=30, dest='maxCost',help = 'Maximum allowable cost.')
    parser.add_argument('--out', type = str, default='out.csv', dest='out', help = 'Output filename')

    args = parser.parse_args()

    #np_directory, np_filename = os.path.split(args.input_npy_filename)
    newframe = buildConnections(args.i1,
     args.i2, args.timepoint1,args.timepoint2,
     greedy=args.greed,openingCost=args.openCost)

    newframe = newframe[newframe.cost < args.maxCost]

    newframe.to_csv(args.out, index = True, header=True)


    return

if __name__ == '__main__':
    main()
