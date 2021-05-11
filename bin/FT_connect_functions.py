# FT_connect_functions

from __future__ import print_function
import numpy as np
import time, os, sys
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage import color, feature, filters, io, measure, morphology, segmentation, img_as_ubyte, transform, registration
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
from shutil import copyfile
from scipy.spatial import distance
from FT_connect_config import *


def genDisplacement(filename_t0,filename_t1):
    global pathToSegs
    temp1 = np.asarray(np.load(filename_t0,allow_pickle=True)).item()
    imfilename1 = temp1['filename'].split('/')[-1]
    img1 = io.imread(pathToSegs+imfilename1);

    temp2 = np.asarray(np.load(filename_t1,allow_pickle=True)).item()
    imfilename2 = temp2['filename'].split('/')[-1]
    img2 = io.imread(pathToSegs+imfilename2);

    shift_vector = registration.phase_cross_correlation(img1, img2)
    return(shift_vector[0])

def buildFeatureFrame(filename_t0,timepoint,pathtoimage="./"):
    temp = np.asarray(np.load(filename_t0,allow_pickle=True)).item()
    imfilename = temp['filename'].split('/')[-1]
    img = io.imread(pathtoimage+imfilename);
    masks = clear_border(temp['masks'])
    image_props = measure.regionprops_table(masks,
                                        intensity_image=img,
                                        properties=('label','area', 'centroid', 'bbox','mean_intensity'))
    im_df = pd.DataFrame(image_props)
    im_df['time'] = timepoint
    return(im_df)

def expandBoundBox(FeatureFrame, expansion = 2):
    hf_row = np.ceil((FeatureFrame['bbox-3']-FeatureFrame['bbox-1'])/2)
    hf_col = np.ceil((FeatureFrame['bbox-2']-FeatureFrame['bbox-0'])/2)
    maxes = expansion*np.amax(np.vstack((hf_row,hf_col)).T,axis=1).astype(int)
    FeatureFrame['ebox-0'] = np.amax(np.vstack((np.zeros(FeatureFrame.shape[0]),FeatureFrame['bbox-0']-maxes)).T,axis=1).astype(int)
    FeatureFrame['ebox-1'] = np.amax(np.vstack((np.zeros(FeatureFrame.shape[0]),FeatureFrame['bbox-1']-maxes)).T,axis=1).astype(int)
    FeatureFrame['ebox-2'] = np.amin(np.vstack((np.zeros(FeatureFrame.shape[0])+np.max(FeatureFrame['bbox-2']),FeatureFrame['bbox-2']+maxes)).T,axis=1).astype(int)
    FeatureFrame['ebox-3'] = np.amin(np.vstack((np.zeros(FeatureFrame.shape[0])+np.max(FeatureFrame['bbox-3']),FeatureFrame['bbox-3']+maxes)).T,axis=1).astype(int)
    return(FeatureFrame)

def futureBoundBox(FeatureFrame,shiftVector):
    FeatureFrame['fbox-0'] = np.amax(np.vstack((np.zeros(FeatureFrame.shape[0]),FeatureFrame['ebox-0'] - shiftVector[1] )).T,axis=1).astype(int)
    FeatureFrame['fbox-1'] = np.amax(np.vstack((np.zeros(FeatureFrame.shape[0]),FeatureFrame['ebox-1'] - shiftVector[0] )).T,axis=1).astype(int)
    FeatureFrame['fbox-2'] = np.amin(np.vstack((np.zeros(FeatureFrame.shape[0])+np.max(FeatureFrame['ebox-2']),FeatureFrame['ebox-2'] - shiftVector[1] )).T,axis=1).astype(int)
    FeatureFrame['fbox-3'] = np.amin(np.vstack((np.zeros(FeatureFrame.shape[0])+np.max(FeatureFrame['ebox-3']),FeatureFrame['ebox-3'] - shiftVector[0] )).T,axis=1).astype(int)
    return(FeatureFrame)

def expectedLocation(FeatureFrame,shiftVector):
    FeatureFrame['fcentroid-0'] = FeatureFrame['centroid-0'] - shiftVector[1]
    FeatureFrame['fcentroid-1'] = FeatureFrame['centroid-1'] - shiftVector[0]
    return(FeatureFrame)

def genCandidateNodes(FeatureFrame_t0, FeatureFrame_t1):
    candidates = (((np.asarray(FeatureFrame_t1['centroid-0'])[:,None]>=np.asarray(FeatureFrame_t0['fbox-0']))&
                   (np.asarray(FeatureFrame_t1['centroid-0'])[:,None]<np.asarray(FeatureFrame_t0['fbox-2']))&
                   (np.asarray(FeatureFrame_t1['centroid-1'])[:,None]>=np.asarray(FeatureFrame_t0['fbox-1']))&
                   (np.asarray(FeatureFrame_t1['centroid-1'])[:,None]<np.asarray(FeatureFrame_t0['fbox-3']))))
    return(candidates)

def getDifference(FeatureFrame_t0, FeatureFrame_t1, feature="position",normed=True):
    if (feature == "position"):
        delta0 = (np.asarray(FeatureFrame_t1['centroid-0'])[:,None]-np.asarray(FeatureFrame_t0['fcentroid-0']))
        delta1 = (np.asarray(FeatureFrame_t1['centroid-1'])[:,None]-np.asarray(FeatureFrame_t0['fcentroid-1']))
        result = np.sqrt(delta0**2 + delta1**2 )
    else :
        result = np.abs(np.asarray(FeatureFrame_t1[feature])[:,None]-np.asarray(FeatureFrame_t0[feature]))
        if normed:
            result = result/10**np.floor(np.log10(np.max(result)))
    return(result)


def DivSizeScore(mom_area, sis_area_1, sis_area_2):
    global DivSizeDiffThreshold,DivSizeScoreReturn,DivSizeRatio_Min,DivSizeRatio_Max

    areaRatio = (mom_area / (sis_area_1 + sis_area_2))
    diffArea = np.abs(sis_area_1 - sis_area_2)/np.sqrt(sis_area_1 *sis_area_2)
    if((areaRatio >= DivSizeRatio_Min)&(areaRatio < DivSizeRatio_Max)&(diffArea<DivSizeDiffThreshold)):
        return(0.0)
    else:
        return(DivSizeScoreReturn)

def DivIntScore(sis_int_1, sis_int_2):
    global DivIntensityDiffThreshold,DivIntensityScoreReturn

    diffInt = np.abs(sis_int_1 - sis_int_2)/np.sqrt(sis_int_1*sis_int_2)
    if((diffInt<DivIntensityDiffThreshold)):
        return(0.0)
    else:
        return(DivIntensityScoreReturn)

def DivScore(FeatureFrame_t0, FeatureFrame_t1, index_mom, index_sis_1, index_sis_2):
    global DivMoveScoreReturn, mitosis_RangeMultiplier

    momFF_select = FeatureFrame_t0.loc[index_mom]
    sis1FF_select = FeatureFrame_t1.loc[index_sis_1]
    sis2FF_select = FeatureFrame_t1.loc[index_sis_2]

    mom_loc = [momFF_select['centroid-0'],momFF_select['centroid-1']]
    mom_corner_1 = [momFF_select['bbox-0'],momFF_select['bbox-1']]
    mom_corner_2 = [momFF_select['bbox-2'],momFF_select['bbox-3']]
    mom_area = (momFF_select['area'])
    mom_range = distance.euclidean(mom_corner_1,mom_corner_2)

    sis1_loc = [sis1FF_select['centroid-0'],sis1FF_select['centroid-1']]
    sis1_area = (sis1FF_select['area'])
    sis1_int = (sis1FF_select['mean_intensity'])

    sis2_loc = [sis2FF_select['centroid-0'],sis2FF_select['centroid-1']]
    sis2_area = (sis2FF_select['area'])
    sis2_int = (sis2FF_select['mean_intensity'])

    mom_s1_dist = distance.euclidean(sis1_loc,mom_loc)
    mom_s2_dist = distance.euclidean(sis2_loc,mom_loc)
    sis_middle_loc = (np.array(sis1_loc)+np.array(sis2_loc))/2

    cost1 = distance.euclidean(sis_middle_loc,mom_loc)
    cost2 = np.abs(mom_s1_dist-mom_s2_dist)
    cost3 = distance.euclidean(sis1_loc,sis2_loc)

    if(cost3 < (mitosis_RangeMultiplier*mom_range)):
        MoveCost = cost1 + cost2/2
    else:
        MoveCost = DivMoveScoreReturn

    SizeCost = DivSizeScore(mom_area=mom_area, sis_area_1=sis1_area, sis_area_2=sis2_area)

    IntCost = DivIntScore(sis_int_1=sis1_int, sis_int_2=sis2_int)

    finalScore = np.round((MoveCost+SizeCost+IntCost),1)
    return([index_mom,index_sis_1,index_sis_2,finalScore])


def GenMitosisPairs(CandidateFrame, motherIndex):
    #returns array of daughter index-pairs in candidate frame
    DaughtersPossible = np.where(CandidateFrame[:,motherIndex])[0]
    if(len(DaughtersPossible)>1):
        DaughtersPairs = np.array(np.meshgrid(DaughtersPossible, DaughtersPossible)).T.reshape(-1,2)
        Sisters = np.unique(np.sort(DaughtersPairs),axis=0)
        Sisters = Sisters[Sisters[:,0] != Sisters[:,1]]
        includeMotherIndex = np.append((np.zeros((Sisters.shape[0],1))+motherIndex).astype(int),Sisters, 1)
        return(includeMotherIndex)
    else:
        return(np.array([[0,0,0]]))

def getMitosisCandidates(CandidateFrame,FeatureFrame_t0, FeatureFrame_t1):
    global mitosis_MaxScore

    divCandidates = np.vstack(list(map(lambda x: GenMitosisPairs(CandidateFrame, x), range(CandidateFrame.shape[1]))))
    divCandidates = divCandidates[(divCandidates[:,1]!=0)&(divCandidates[:,2]!=0)]
    divScores = np.vstack(list(map(lambda x: DivScore(FeatureFrame_t0=FeatureFrame_t0, FeatureFrame_t1=FeatureFrame_t1,
                                                      index_mom=divCandidates[x,0], index_sis_1=divCandidates[x,1],
                                                      index_sis_2=divCandidates[x,2]),
                                   range(divCandidates.shape[0]))))
    divScores = divScores[divScores[:,3]<mitosis_MaxScore]
    return(divScores[:,:3].astype(int))

def getCostMatrix(FeatureFrame_t0, FeatureFrame_t1, shiftVec):
    global track_frameExpCoeff, costIntCoefficient, costSizeCoefficient, costPositionCoefficient

    FeatureFrame_t0 = expandBoundBox(FeatureFrame_t0, expansion = track_frameExpCoeff)
    FeatureFrame_t0 = futureBoundBox(FeatureFrame_t0, shiftVector=shiftVec)
    CandidateMtx = genCandidateNodes(FeatureFrame_t0,FeatureFrame_t1)
    FeatureFrame_t0 = expectedLocation(FeatureFrame_t0, shiftVector=shiftVec)
    deltaPosition = getDifference(FeatureFrame_t0,FeatureFrame_t1,"position")
    deltaArea = getDifference(FeatureFrame_t0,FeatureFrame_t1,"area")
    deltaIntensity = getDifference(FeatureFrame_t0,FeatureFrame_t1,"mean_intensity")

    costMatrix = ((costIntCoefficient*(deltaIntensity)+costSizeCoefficient*(deltaArea)+costPositionCoefficient*(deltaPosition))*CandidateMtx)
    return((FeatureFrame_t0,CandidateMtx,costMatrix))

def solveMinCostFlow(CostMatrix,mitosisCands):
    global openingCost, closingCost

    t0Nodes = np.array(range(CostMatrix.shape[1]))+1
    t1Nodes = np.array(range(CostMatrix.shape[1],np.sum(CostMatrix.shape)))+1

    start_nodes = np.concatenate((np.repeat([0],np.sum(CostMatrix.shape)), # all connections to source node
                                  (np.nonzero(CostMatrix))[1]+1, # all connections from t0 to t1
                                  t0Nodes,t1Nodes # all connections to sink node
                                 )).tolist()

    end_nodes = np.concatenate((t0Nodes,t1Nodes, # all connections to source node
                                np.nonzero(CostMatrix)[0]+int(CostMatrix.shape[1])+1, # all connections from t0 to t1
                                np.repeat([np.sum(CostMatrix.shape)+1],np.sum(CostMatrix.shape)) # all connections to sink node
                                 )).tolist()

    costs = np.concatenate((np.repeat(1,CostMatrix.shape[1]), # all connections to source node
                            np.repeat(openingCost,CostMatrix.shape[0]),
                            CostMatrix[CostMatrix!=0], # all connections from t0 to t1
                            np.repeat(closingCost,CostMatrix.shape[1]), # all connections to sink node
                            np.repeat(1,CostMatrix.shape[0]) # all connections to sink node
                            )).tolist()

    nodeCaps = np.concatenate((t0Nodes,t1Nodes),axis=0)
    nodeCaps = np.vstack((nodeCaps, np.repeat(1,len(nodeCaps)))).T
    if(len(mitosisCands)>0):
        nodeCaps[np.searchsorted(nodeCaps[:,0],mitosisCands+1),1]=2

    capacities = np.concatenate((nodeCaps[:,1], # all connections to source node
                                 np.repeat(1,np.sum(CostMatrix[CostMatrix!=0].shape)),# all connections from t0 to t1
                                np.repeat(1,np.sum(CostMatrix.shape)) # all connections to sink node
                                )).tolist()

#    supply_amount = np.min([CostMatrix.shape[1]+len(mitosisCands),CostMatrix.shape[0]])#np.max([CostMatrix.shape[0],CostMatrix.shape[1]])
    supply_amount = np.max([CostMatrix.shape[1],CostMatrix.shape[0]])#np.max([CostMatrix.shape[0],CostMatrix.shape[1]])
    supplies = np.concatenate(([supply_amount],np.repeat(0,np.sum(CostMatrix.shape)),[-1*supply_amount])).tolist()


    min_cost_flow = pywrapgraph.SimpleMinCostFlow()

    # Add each arc.
    for i in range(len(start_nodes)):
        min_cost_flow.AddArcWithCapacityAndUnitCost(start_nodes[i], end_nodes[i],capacities[i], int(costs[i]))

    for i in range(len(supplies)):
        min_cost_flow.SetNodeSupply(i, supplies[i])


    ArcFrame = pd.DataFrame()
    # Find the minimum cost flow between node 0 and node 4.
    if min_cost_flow.Solve() == min_cost_flow.OPTIMAL:
        print('Minimum cost:', min_cost_flow.OptimalCost())
        for i in range(min_cost_flow.NumArcs()):
            cost = min_cost_flow.Flow(i) * min_cost_flow.UnitCost(i)

            ArcFrame = ArcFrame.append(pd.DataFrame([min_cost_flow.Tail(i),
              min_cost_flow.Head(i),
              min_cost_flow.Flow(i),
              min_cost_flow.Capacity(i),
              cost]).T)
    else:
        print('There was an issue with the min cost flow input.')

    ArcFrame = ArcFrame.rename(columns={0:'start',1:'end',2:"Flow",3:"Capacity",4:"Cost"})
    FinalFrame = ArcFrame.loc[ArcFrame["Flow"]!=0,]
    FinalFrame = FinalFrame[(FinalFrame["start"]!=0)&(FinalFrame["end"]!=(np.sum(CostMatrix.shape)+1))]
    FinalFrame['end']=(FinalFrame['end']-int(CostMatrix.shape[1])-1)
    FinalFrame['start']=(FinalFrame['start']-1)
    FinalFrame = FinalFrame.reset_index(drop=True)
    return(FinalFrame)

def formatFinalConnections(Connections,FeatureFrame_t0,FeatureFrame_t1,timept_0, timept_1):

    t0_connected_labels = np.array(FeatureFrame_t0.loc[Connections['start']]['label'])
    t1_connected_labels = np.array(FeatureFrame_t1.loc[Connections['end']]['label'])
    t0_die_lables = np.array(FeatureFrame_t0[~FeatureFrame_t0.index.isin(np.array(Connections['start']))]['label'])
    t1_born_lables = np.array(FeatureFrame_t1[~FeatureFrame_t1.index.isin(np.array(Connections['end']))]['label'])

    connectedFrames = pd.DataFrame([t0_connected_labels,t1_connected_labels]).T
    birthFrames = pd.DataFrame([np.repeat(-1,len(t1_born_lables)),t1_born_lables]).T
    deathFrames = pd.DataFrame([t0_die_lables,np.repeat(-1,len(t0_die_lables))]).T
    connectedFrames = connectedFrames.rename(columns={0: "start", 1: "end"})
    birthFrames = birthFrames.rename(columns={0: "start", 1: "end"})
    deathFrames = deathFrames.rename(columns={0: "start", 1: "end"})

    u, c = np.unique(connectedFrames['start'], return_counts=True)
    dup = u[c >1]
    splitFrames = connectedFrames[connectedFrames.start.isin(dup)]
    nonsplitFrames = connectedFrames[~connectedFrames.start.isin(dup)]
    if (splitFrames.shape[0]>0):
        splitFrames.loc[:,"annotation"]="split"
    if (birthFrames.shape[0]>0):
        birthFrames.loc[:,"annotation"]="birth"
    if (deathFrames.shape[0]>0):
        deathFrames.loc[:,"annotation"]="death"
    if (nonsplitFrames.shape[0]>0):
        nonsplitFrames.loc[:,"annotation"]="pass"

    DetailFrames = pd.concat([nonsplitFrames,splitFrames,birthFrames,deathFrames])
    DetailFrames = DetailFrames.reset_index(drop=True)
    DetailFrames['t_start']=timept_0
    DetailFrames['t_end']=timept_1
    DetailFrames['MasterID_'+timept_0]=[str(x) + '_' + str(y) for x, y in zip(DetailFrames['t_start'], DetailFrames['start'])]
    DetailFrames['MasterID_'+timept_1]=[str(x) + '_' + str(y) for x, y in zip(DetailFrames['t_end'], DetailFrames['end'])]
    DetailFrames['rlgID_'+timept_0]=[(10**5)*(int(x)+1) + int(y) for x, y in zip(DetailFrames['t_start'], DetailFrames['start'])]
    DetailFrames['rlgID_'+timept_1]=[(10**5)*(int(x)+1) + int(y) for x, y in zip(DetailFrames['t_end'], DetailFrames['end'])]

    return(DetailFrames)

def RLGwrap(finalConnections,FeatureFrame_t0,FeatureFrame_t1,time0=0):
    ffC_nodeath = finalConnections[finalConnections['annotation']!='death']
    ffC_nodeath['xCoord']=np.array((FeatureFrame_t1.iloc[pd.Index(FeatureFrame_t1['label']).get_indexer(ffC_nodeath['end'])]['centroid-1']).astype(int))
    ffC_nodeath['yCoord']=np.array((FeatureFrame_t1.iloc[pd.Index(FeatureFrame_t1['label']).get_indexer(ffC_nodeath['end'])]['centroid-0']).astype(int))
    ffC_nodeath =ffC_nodeath.rename(columns={'rlgID_'+str(int(np.unique(ffC_nodeath['t_end']))):'pointID',
                                             'rlgID_'+str(int(np.unique(ffC_nodeath['t_start']))):'parentID'})
    ffC_nodeath['frameNumber']=time0+2

    ffC_rlg = ffC_nodeath[['pointID','frameNumber','xCoord','yCoord','parentID','annotation']]
    ffC_rlg.loc[:,'parentID'] = ffC_rlg['parentID'].astype(str)
    ffC_rlg.loc[ffC_rlg['annotation']=="birth",'parentID'] = 'NaN'

    if(time0==0):
        ffC_nobirth = finalConnections[(finalConnections['annotation']!='birth')]
        ffC_nobirth['pointID']=ffC_nobirth['rlgID_'+str(time0)]
        ffC_nobirth['parentID'] = 'NaN'

        ffC_nobirth['xCoord']=np.array((FeatureFrame_t0.iloc[pd.Index(FeatureFrame_t0['label']).get_indexer(ffC_nobirth['start'])]['centroid-1']).astype(int))
        ffC_nobirth['yCoord']=np.array((FeatureFrame_t0.iloc[pd.Index(FeatureFrame_t0['label']).get_indexer(ffC_nobirth['start'])]['centroid-0']).astype(int))
        ffC_nobirth['frameNumber']=(time0+1)

        ffC_nobirth = ffC_nobirth[['pointID','frameNumber','xCoord','yCoord','parentID','annotation']].drop_duplicates()

        ffC_rlg = pd.concat([ffC_nobirth,ffC_rlg])

    return(ffC_rlg)

def fullPipeline(file_1,timept_1,file_2,timept_2):
    global pathToSegs,rlg_out
    global track_frameExpCoeff,track_openingCost,track_closingCost,mitosis_RangeMultiplier,mitosis_MaxScore

    shifts = genDisplacement(pathToSegs+file_1,pathToSegs+file_2)
    t0 = buildFeatureFrame(pathToSegs+file_1,timept_1,pathToSegs)
    t1 = buildFeatureFrame(pathToSegs+file_2,timept_2,pathToSegs)

    t0, candies, costMatrix1 = getCostMatrix(t0, t1,shiftVec =shifts)

    mitosisCandidates =getMitosisCandidates(CandidateFrame=candies,FeatureFrame_t0=t0,FeatureFrame_t1=t1)

    if (mitosisCandidates.shape[0]>1):
        mitosisConnections = np.unique(np.append(np.delete(mitosisCandidates, 1, 1),
                                                 np.delete(mitosisCandidates, 2, 1),axis=0),axis=0)

        for i in range(mitosisConnections.shape[0]):
            costMatrix1[mitosisConnections[i,1],mitosisConnections[i,0]] = (costMatrix1[mitosisConnections[i,1],mitosisConnections[i,0]]/2)

        mitosisCandidateNodes = np.unique(mitosisConnections[:,0])

        FinFrame = solveMinCostFlow(CostMatrix=costMatrix1, mitosisCands=mitosisCandidateNodes)
    else:
        FinFrame = solveMinCostFlow(CostMatrix=costMatrix1,mitosisCands=[])


    final = formatFinalConnections(Connections=FinFrame,FeatureFrame_t0=t0,FeatureFrame_t1=t1,timept_0=timept_1, timept_1=timept_2)
    return(t0,t1,final)
