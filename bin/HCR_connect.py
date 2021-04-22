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

def buildFeatureFrame(filename,timepoint):
    temp = np.asarray(np.load(filename,allow_pickle=True)).item()
    imfilename = temp['filename']
    img = io.imread(imfilename);
    masks = clear_border(temp['masks'])
    image_props = measure.regionprops_table(masks,
                                        intensity_image=img,
                                        properties=('label','area','filled_area', 'centroid',
                                                    'eccentricity','mean_intensity'))
    im_df = pd.DataFrame(image_props)
    im_df['time'] = timepoint
    return(im_df)

def generateCandidates(image1, image2, im1_select, dist_multiplier=2):
    delX = np.sqrt((image1['centroid-0'][im1_select]-image2['centroid-0'])**2+
                   (image1['centroid-1'][im1_select]-image2['centroid-1'])**2)
    max_dist = dist_multiplier*min(delX)
    candidates = np.array(delX[delX < max_dist].index)
    return(candidates)


def generateLinks(filename_t0, filename_t1,timepoint, nnDist = 10,costMax=35, mN_Int = 10, mN_Ecc=4, mN_Area=25, mN_Disp=1):
    ip0 = buildFeatureFrame(filename_t0,timepoint)
    ip1 = buildFeatureFrame(filename_t1,timepoint+1)

    arr = pd.DataFrame()

    for i in np.array(ip0.index):
        candidates = generateCandidates(ip0, ip1, i, dist_multiplier=nnDist)
        canFRAME = pd.DataFrame(candidates)
        canFRAME["1"] = i
        arr = arr.append(canFRAME)

    arr = arr.rename(columns={0: "t1", "1": "t0"})
    arr = arr.reset_index(drop=True)

    properties = pd.DataFrame()

    mInt_0 = float(np.median(ip0.loc[:,['mean_intensity']]))
    mInt_1 = float(np.median(ip1.loc[:,['mean_intensity']]))

    for link in np.array(arr.index):
        tmp_props_0 = (ip0.loc[arr.loc[link,["t0"]],:])
        tmp_props_1 = (ip1.loc[arr.loc[link,["t1"]],:])
        deltaInt = (np.abs((int(tmp_props_0["mean_intensity"])/mInt_0)-(int(tmp_props_1["mean_intensity"])/mInt_1))/
                np.mean([(int(tmp_props_0["mean_intensity"])/mInt_0),(int(tmp_props_1["mean_intensity"])/mInt_1)]))
        deltaArea = (np.abs(int(tmp_props_0['area']) - int(tmp_props_1['area']))/
                np.mean([int(tmp_props_0["area"]),int(tmp_props_1["area"])]))
        deltaEcc = np.absolute(float(tmp_props_0['eccentricity']) - float(tmp_props_1['eccentricity']))
        deltaX = np.sqrt((int(tmp_props_0['centroid-0'])-int(tmp_props_1['centroid-0']))**2+
                     (int(tmp_props_0['centroid-1'])-int(tmp_props_1['centroid-1']))**2)
        properties = properties.append(pd.DataFrame([int(tmp_props_0['label']),int(tmp_props_1['label']),
                            deltaInt ,deltaArea,deltaEcc,deltaX]).T)

    properties = properties.rename(columns={0: "label_t0", 1: "label_t1", 2: "deltaInt",
                                           3: "deltaArea", 4: "deltaEcc", 5: "deltaX"})
    properties = properties.reset_index(drop=True)
    properties["Cost"]=(properties.loc[:,"deltaInt"]*mN_Int)+(properties.loc[:,"deltaEcc"]*mN_Ecc)+(properties.loc[:,"deltaArea"]*mN_Area)+(properties.loc[:,"deltaX"]*mN_Disp)
    properties["TransitionCapacity"]=1
    properties = properties.loc[properties["Cost"]<costMax]
    properties = properties.reset_index(drop=True)

    return(properties)

def DivSimScore(daughterCell_1, daughterCell_2, FrameNext):
    daughterStats_1 = FrameNext[(FrameNext['label'] == daughterCell_1)]
    daughterStats_2 = FrameNext[(FrameNext['label'] == daughterCell_2)]

    deltaInt = (np.abs((int(daughterStats_1["mean_intensity"]))-(int(daughterStats_2["mean_intensity"])))/
            np.mean([(int(daughterStats_1["mean_intensity"])),(int(daughterStats_2["mean_intensity"]))]))

    deltaArea = (np.abs(int(daughterStats_1['area']) - int(daughterStats_2['area']))/
            np.mean([int(daughterStats_1["area"]),int(daughterStats_2["area"])]))

    deltaEcc = np.absolute(float(daughterStats_1['eccentricity']) - float(daughterStats_2['eccentricity']))

    deltaX = np.sqrt((int(daughterStats_1['centroid-0'])-int(daughterStats_2['centroid-0']))**2+
                     (int(daughterStats_1['centroid-1'])-int(daughterStats_2['centroid-1']))**2)

    sims = pd.DataFrame([int(daughterCell_1),int(daughterCell_2),
                         deltaInt ,deltaArea,deltaEcc,deltaX]).T
    sims = sims.rename(columns={0: "label_D1", 1: "label_D2", 2: "D2deltaInt",
                                           3: "D2deltaArea", 4: "D2deltaEcc", 5: "D2deltaX"})
    return(sims)

def DivSetupScore(motherCell, daughterCell_1, daughterCell_2, FrameCurr, FrameNext):
    #determine similarities between mother and daughters
    simDF = DivSimScore(daughterCell_1, daughterCell_2, FrameNext)

    #determine relative area of mother compared to daughters
    MotherArea = int(FrameCurr[(FrameCurr['label'] == motherCell)]['area'])
    daughterArea_1 = int(FrameNext[(FrameNext['label'] == daughterCell_1)]['area'])
    daughterArea_2 = int(FrameNext[(FrameNext['label'] == daughterCell_2)]['area'])
    areaChange = MotherArea/(daughterArea_1 + daughterArea_2)

    simDF["MDDeltaArea"] = areaChange

    return(simDF)

def DivisionCanditates(propMtx, filename_t0,filename_t1,timepoint,mS_Area = 10, mS_Ecc = 2, mS_Int = 2, mS_Disp = 1, MDAR_thresh = 0.75, SDis_thresh = 20.0):
    ip0 = buildFeatureFrame(filename_t0,timepoint)
    ip1 = buildFeatureFrame(filename_t1,timepoint+1)

    Mothers = np.unique(propMtx.loc[:,['label_t0']])
    DivCandidacy = pd.DataFrame()

    for cell in Mothers:
        DaughtersPossible = (propMtx[(propMtx['label_t0'] == cell)].loc[:,'label_t1'])
        DaughtersPairs = np.array(np.meshgrid(DaughtersPossible, DaughtersPossible)).T.reshape(-1,2)
        Sisters = np.unique(np.sort(DaughtersPairs),axis=0)
        for pair in range(Sisters.shape[0]):
            if (Sisters[pair,0] != Sisters[pair,1]):
                tmpScoreSetup = (DivSetupScore(cell,Sisters[pair,0], Sisters[pair,1], ip0,ip1))
                LogicMDAR = (tmpScoreSetup["MDDeltaArea"]>MDAR_thresh)
                ScoreSDis = (mS_Int*tmpScoreSetup["D2deltaInt"]) + (mS_Area*tmpScoreSetup["D2deltaArea"]) + (mS_Ecc*tmpScoreSetup["D2deltaEcc"]) + (mS_Disp*tmpScoreSetup["D2deltaX"])
                LogicSDis = (ScoreSDis<SDis_thresh)
                tmpCandidacy = pd.DataFrame([cell,Sisters[pair,0],Sisters[pair,1],(LogicSDis&LogicMDAR).bool()]).T
                DivCandidacy = DivCandidacy.append(tmpCandidacy)

    DivCandidacy = DivCandidacy.rename(columns={0: "Mother", 1: "Daughter1", 2: "Daughter2",3: "Div"})
    DivCandidacy = DivCandidacy.reset_index(drop=True)

    # select true values
    DivSelect = DivCandidacy[(DivCandidacy['Div'] == True)]
    DivConnects_1 = DivSelect[['Mother','Daughter1','Div']]
    DivConnects_2 = DivSelect[['Mother','Daughter2','Div']]
    DivConnects_1 = DivConnects_1.rename(columns={'Mother': "label_t0", 'Daughter1': "label_t1"})
    DivConnects_2 = DivConnects_2.rename(columns={'Mother': "label_t0", 'Daughter2': "label_t1"})
    DivConnects = pd.concat([DivConnects_1,DivConnects_2])
    DivConnects = DivConnects.reset_index(drop=True)

    return(DivConnects)

def UpdateConnectionsDiv(propMtx,DivCandidatesMtx):
    propMtx.loc[propMtx['label_t0'].isin(np.unique(DivCandidatesMtx['label_t0'])),['TransitionCapacity']] = 2

    for div in range(DivCandidatesMtx.shape[0]):
        tmp_prop = propMtx.loc[(DivCandidatesMtx.loc[div,'label_t0'] ==propMtx['label_t0'])&(DivCandidatesMtx.loc[div,'label_t1'] ==propMtx['label_t1']),]
        old_score = float(tmp_prop.loc[:,'Cost'])
        new_score = (old_score/2)
        propMtx.loc[(DivCandidatesMtx.loc[div,'label_t0'] ==propMtx['label_t0'])&(DivCandidatesMtx.loc[div,'label_t1'] ==propMtx['label_t1']),'Cost'] = new_score

    return(propMtx)

def SolveMinCostTable(filename_t0, filename_t1, DivisionTable,timepoint, OpeningCost = 30, ClosingCost = 30):
    #rename
    ip0 = buildFeatureFrame(filename_t0,timepoint)
    ip0 = ip0.rename(columns={"label" : "label_t0"})
    ip1 = buildFeatureFrame(filename_t1,timepoint+1)
    ip1 = ip1.rename(columns={"label" : "label_t1"})

    ip0["slabel_t0"] = np.array(range(ip0.label_t0.shape[0]))+1
    i0max = np.max(np.asarray(ip0["slabel_t0"]))
    ip1["slabel_t1"] = np.array(range(i0max,i0max+ip1.label_t1.shape[0]))+1
    i1max = np.max(np.asarray(ip1["slabel_t1"]))

    i0_translation = ip0[["label_t0","slabel_t0"]]
    i1_translation = ip1[["label_t1","slabel_t1"]]

    result_tmp = pd.merge(DivisionTable, i0_translation, on=['label_t0'])
    result = pd.merge(result_tmp, i1_translation, on=['label_t1'])
    result_shorthand = result[['slabel_t0','slabel_t1','Cost','TransitionCapacity']]

    transNodes0 = np.array(result_shorthand['slabel_t0']) ;
    transNodes1 = np.array(result_shorthand['slabel_t1']) ;
    transCosts = np.array(result_shorthand['Cost']) ;
    transCaps = np.repeat(1,transNodes0.size) ;


    sourceNodes0 = np.repeat([0],i1max)
    sourceNodes1 = np.array(range(i1max))+1
    sourceCosts = np.concatenate((np.repeat(1,ip0.shape[0]),np.repeat(OpeningCost,ip1.shape[0])), axis=None)
    #Source capacities are dictates by which node could be splitting. Source capacity = 2 if there was a division candidate
    tmpUnique0 = result_shorthand[["slabel_t0","TransitionCapacity"]].drop_duplicates()
    HighCaps = tmpUnique0.loc[tmpUnique0["TransitionCapacity"]==2,]
    LowCaps = pd.DataFrame(i0_translation).copy(deep=True)
    LowCaps['Cap'] = 1
    LowCaps.loc[LowCaps['slabel_t0'].isin(np.array(HighCaps['slabel_t0'])),'Cap'] = 2
    sourceCaps = np.concatenate((np.array(LowCaps['Cap']),np.repeat(1,ip1.shape[0])), axis=None)


    sinkNodes0 = np.array(range(i1max))+1
    sinkNodes1 = np.repeat([i1max+1],i1max)
    sinkCosts = np.concatenate((np.repeat(ClosingCost,ip0.shape[0]),np.repeat(1,ip1.shape[0])), axis=None)
    sinkCaps = np.repeat(1,i1max)

    # Define the directed graph for the flow.
    min_cost_flow = pywrapgraph.SimpleMinCostFlow()

    start_nodes = np.concatenate((sourceNodes0, transNodes0, sinkNodes0)).tolist()
    end_nodes = np.concatenate((sourceNodes1, transNodes1, sinkNodes1)).tolist()
    capacities = np.concatenate((sourceCaps, transCaps, sinkCaps)).tolist()
    costs  = np.concatenate((sourceCosts, transCosts, sinkCosts)).tolist()
    source = 0
    sink = i1max+1
    supply_amount = np.max([i0max,i1max-i0max])
    supplies = np.concatenate(([supply_amount],np.repeat(0,i1max),[-1*supply_amount])).tolist()

    min_cost_flow = pywrapgraph.SimpleMinCostFlow()

    # Add each arc.
    for i in range(len(start_nodes)):
        min_cost_flow.AddArcWithCapacityAndUnitCost(start_nodes[i], end_nodes[i],capacities[i], int(costs[i]))
      # Add node supplies.

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
    #ArcFrame = ArcFrame.reset_index(drop=True)
    FinalFrame = ArcFrame.loc[ArcFrame["Flow"]!=0,]
    FinalFrame = FinalFrame.reset_index(drop=True)

    return(FinalFrame)

def ReviewCostTable(minCostFlowtable, timepoint, OpeningCost=30,ClosingCost=30):
    sink = max(minCostFlowtable["end"])
    Transitions = minCostFlowtable.loc[(minCostFlowtable["start"]!=0)&(minCostFlowtable["end"]!=sink),]

    trans_start_nodes = np.unique(Transitions["start"])
    trans_end_nodes = np.unique(Transitions["end"])

    #find nodes that either appear (no start) or disappear (no end)
    appearing = minCostFlowtable[(~minCostFlowtable.start.isin(trans_start_nodes))&
                                 (~minCostFlowtable.end.isin(trans_start_nodes))&
                                 (~minCostFlowtable.start.isin(trans_end_nodes))&
                                 (~minCostFlowtable.end.isin(trans_end_nodes))]
    appearing = appearing.loc[(appearing["Cost"] == OpeningCost)|(appearing["Cost"] == ClosingCost)]
    appearing = appearing.reset_index(drop=True)
    appearFrame = pd.DataFrame()
    for i in range(appearing.shape[0]):
        if(appearing.loc[i,"start"] == 0):
            appearFrame = appearFrame.append(pd.DataFrame([-1,appearing.loc[i,"end"]]).T)
        elif(appearing.loc[i,"end"] == sink):
            appearFrame = appearFrame.append(pd.DataFrame([appearing.loc[i,"end"],-1]).T)
    appearFrame = appearFrame.rename(columns={0:"slabel_t0",1:"slabel_t1"})
    appearFrame = appearFrame.reset_index(drop=True)

    #Assemble
    transFrame = Transitions.loc[:,["start","end"]]
    transFrame = transFrame.rename(columns={"start":"slabel_t0","end":"slabel_t1"})
    totalFrame = pd.concat([appearFrame,transFrame])
    totalFrame = totalFrame.reset_index(drop=True)
    totalFrame["timepoint"] = timepoint

    return(totalFrame)


def TranslationTable(filename_t0, filename_t1, DivisionTable,timepoint):
    #rename
    ip0 = buildFeatureFrame(filename_t0,timepoint)
    ip0 = ip0.rename(columns={"label" : "label_t0"})
    ip1 = buildFeatureFrame(filename_t1,timepoint+1)
    ip1 = ip1.rename(columns={"label" : "label_t1"})

    ip0["slabel_t0"] = np.array(range(ip0.label_t0.shape[0]))+1
    i0max = np.max(np.asarray(ip0["slabel_t0"]))
    ip1["slabel_t1"] = np.array(range(i0max,i0max+ip1.label_t1.shape[0]))+1
    i1max = np.max(np.asarray(ip1["slabel_t1"]))

    i0_translation = ip0[["label_t0","slabel_t0"]]
    i1_translation = ip1[["label_t1","slabel_t1"]]

    dvtabDF = DivisionTable
    result_tmp = pd.merge(dvtabDF, i0_translation, on=['label_t0'])
    translation_table = pd.merge(result_tmp, i1_translation, on=['label_t1'])
    #result_shorthand = result[['slabel_t0','slabel_t1','Cost','TransitionCapacity']]

    startLabels = translation_table.loc[:,["label_t0","slabel_t0"]]
    startLabels["timepoint"] = timepoint
    startLabels["frame"] = timepoint+1

    endLabels = translation_table.loc[:,["label_t1","slabel_t1"]]
    endLabels["timepoint"] = timepoint+1
    endLabels["frame"] = timepoint+2

    startLabels = startLabels.rename(columns={"label_t0":"label","slabel_t0":"slabel"})
    endLabels = endLabels.rename(columns={"label_t1":"label","slabel_t1":"slabel"})
    allLabels = pd.concat([startLabels,endLabels])
    allLabels = allLabels.reset_index(drop=True)
    allLabels = allLabels.astype( 'int64')

    allLabels["Master_ID"] = allLabels["timepoint"].astype('str')+"_"+allLabels["label"].astype('str')
    allLabels = allLabels.astype({"Master_ID":'str'})

    allLabels["RajTLG_ID"] = allLabels["frame"]*int(10**(np.ceil(np.log10(max(allLabels['slabel'])))+2))+allLabels["label"]
    allLabels = allLabels.drop_duplicates()
    allLabels = allLabels.reset_index(drop=True)
    return(allLabels)


def TranslateConnections(ConnectionTable, TranslationTable, timepoint, preference = "Master_ID"):

    subTranslationTable_0 = TranslationTable.loc[:,[preference,"slabel"]]
    subTranslationTable_0['slabel_t0'] = subTranslationTable_0['slabel']
    subTranslationTable_1 = TranslationTable.loc[:,[preference,"slabel"]]
    subTranslationTable_1['slabel_t1'] = subTranslationTable_1['slabel']

    merge_0 = pd.merge(ConnectionTable, subTranslationTable_0, on="slabel_t0")
    merge = pd.merge(merge_0, subTranslationTable_1, on="slabel_t1")

    pref = str(preference)

    result = merge.loc[:,[pref+"_x",pref+"_y"]]
    result = result.drop_duplicates()
    result =  result.dropna(thresh=1)
    result = result.reset_index(drop=True)
    result = result.rename(columns = {(pref+"_x") : (pref+"_"+str(timepoint)), (pref+"_y") : (pref+"_"+str(timepoint+1))})
    return(result)

def RajTLG_wrap(filename_t0, filename_t1,timepoint,ConnectionTable,TranslationTable):
    frame0 = buildFeatureFrame(filename_t0,timepoint);
    frame1 = buildFeatureFrame(filename_t1,timepoint+1);
    frames = pd.concat([frame0,frame1])
    frames["timepoint"] = frames["time"]

    InfoDF = pd.merge(frames,TranslationTable, on=['label','timepoint'])

    RajTLG_translation = TranslateConnections(ConnectionTable=ConnectionTable, TranslationTable=TranslationTable, timepoint=timepoint, preference="RajTLG_ID")

    RajTLGFrame = pd.DataFrame()
    if (timepoint == 0):
        for i in range(RajTLG_translation.shape[0]):
            tmpID = RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint)]
            tmpFrame = int(InfoDF.loc[InfoDF["RajTLG_ID"] == RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint)],"frame"])
            tmpX = int(InfoDF.loc[InfoDF["RajTLG_ID"] == RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint)],"centroid-1"])
            tmpY = int(InfoDF.loc[InfoDF["RajTLG_ID"] == RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint)],"centroid-0"])
            tmpParent = "NaN"
            RajTLGFrame = RajTLGFrame.append(pd.DataFrame([tmpID,tmpFrame,tmpX,tmpY,tmpParent]).T)

    for i in range(RajTLG_translation.shape[0]):
            tmpID = RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint+1)]
            tmpFrame = int(InfoDF.loc[InfoDF["RajTLG_ID"] == RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint+1)],"frame"])
            tmpX = int(InfoDF.loc[InfoDF["RajTLG_ID"] == RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint+1)],"centroid-1"])
            tmpY = int(InfoDF.loc[InfoDF["RajTLG_ID"] == RajTLG_translation.loc[i,"RajTLG_ID"+"_"+str(timepoint+1)],"centroid-0"])
            tmpParent = int(RajTLG_translation.loc[RajTLG_translation["RajTLG_ID"+"_"+str(timepoint+1)] == tmpID,
                           "RajTLG_ID"+"_"+str(timepoint)])
            RajTLGFrame = RajTLGFrame.append(pd.DataFrame([tmpID,tmpFrame,tmpX,tmpY,tmpParent]).T)

    RajTLGFrame = RajTLGFrame.reset_index(drop=True)
    RajTLGFrame = RajTLGFrame.rename(columns={0:"pointID", 1:"frameNumber",
                                             2:"xCoord",3:"yCoord",4:"parentID"})
    RajTLGFrame["annotation"] = "none"
    #RajTLGFrame.to_csv(outfilename,index=False)
    return(RajTLGFrame)


def HCR_connect(sampleName, TLlast_mask, HCR_mask, timepoint, nnDist=3, costMax=35, mN_Int=10, mN_Ecc=4, mN_Area=25, mN_Disp=1, mS_Area = 10, mS_Ecc = 2, mS_Int = 2, mS_Disp = 1, MDAR_thresh = 0.75, SDis_thresh = 20.0, openingCost = 30, closingCost = 30):
    propies = generateLinks(filename_t0 = TLlast_mask, filename_t1 = HCR_mask,
                            timepoint = timepoint, nnDist = nnDist,
                            costMax = costMax, mN_Int = mN_Int,
                            mN_Ecc = mN_Ecc, mN_Area = mN_Area,
                            mN_Disp = mN_Disp)

    tmpdivs = DivisionCanditates(propMtx = propies,
                                 filename_t0 = TLlast_mask, filename_t1 = HCR_mask,
                                 MDAR_thresh = MDAR_thresh, SDis_thresh = SDis_thresh,
                                 mS_Disp = mS_Disp, mS_Area = mS_Area,
                                 mS_Ecc = mS_Ecc, mS_Int = mS_Int,
                                 timepoint = timepoint)

    finaldivs = UpdateConnectionsDiv(propies, tmpdivs)

    minCost_table = SolveMinCostTable(TLlast_mask, HCR_mask,
                                      DivisionTable=finaldivs,
                                      timepoint=timepoint,
                                      OpeningCost = openingCost,
                                      ClosingCost = closingCost)

    finTable = ReviewCostTable(minCostFlowtable = minCost_table, timepoint=timepoint)

    translation_table = TranslationTable(TLlast_mask, HCR_mask, DivisionTable=finaldivs,
                                         timepoint=timepoint)

    masterConnects_Raj = TranslateConnections(finTable, translation_table, timepoint=timepoint, preference="RajTLG_ID")
    masterConnects_Master = TranslateConnections(finTable, translation_table, timepoint=timepoint, preference="Master_ID")

    col_df = finTable[(finTable['slabel_t0']!=-1)&(finTable['slabel_t1']!=-1)]

    col_df.to_csv('results/'+sampleName+'/HCR/'+sampleName+'_HCR_connect.csv', index=False)
    translation_table.to_csv('results/'+sampleName+'/HCR/'+sampleName+'_HCR_translation.csv', index=False)
    masterConnects_Raj.to_csv('results/'+sampleName+'/HCR/'+sampleName+'_HCR_connections_RajLab.csv', index=False)
    masterConnects_Master.to_csv('results/'+sampleName+'/HCR/'+sampleName+'_HCR_connections_MasterID.csv', index=False)
