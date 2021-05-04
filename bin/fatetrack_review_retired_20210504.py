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
from fatetrack_connections import buildFeatureFrame, buildOffsetFrame, generateCandidates, generateLinks, DivSimScore, DivSetupScore, DivisionCanditates, UpdateConnectionsDiv, TranslationTable, SolveMinCostTable, ReviewCostTable

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

def RajTLG_wrap(filename_t0, filename_t1,timepoint,ConnectionTable,TranslationTable,path="./"):
    frame0 = buildFeatureFrame(filename_t0,timepoint,pathtoimage=path);
    frame1 = buildFeatureFrame(filename_t1,timepoint+1,pathtoimage=path);
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

def MatchToGoldStd(FileCompare,FileGoldSTD):
    GoldSTD = pd.read_csv(FileGoldSTD)
    FateTrack = pd.read_csv(FileCompare)

    GoldTranslationTable = pd.DataFrame()

    for obj in range(FateTrack.shape[0]):
        FateID = FateTrack.loc[obj,"pointID"]
        frame = FateTrack.loc[obj,"frameNumber"]
        xC = FateTrack.loc[obj,"xCoord"]
        yC = FateTrack.loc[obj,"yCoord"]
        tmpGold = GoldSTD.loc[GoldSTD["frameNumber"] == frame,]
        tmpGold = tmpGold.reset_index(drop=True)
        dist = np.array(np.sqrt((tmpGold["xCoord"]-xC)**2 + (tmpGold["yCoord"]-yC)**2))
        GoldIndex = np.where(dist == dist.min())[0][0]
        GoldID = tmpGold.loc[GoldIndex,"pointID"]
        GoldTranslationTable = GoldTranslationTable.append(pd.DataFrame([GoldID,FateID]).T)

    GoldTranslationTable = GoldTranslationTable.rename(columns={0:"GoldID",1:"FateID"})
    return(GoldTranslationTable)

def CheckAccuracy(frame,FileCompare,FileGoldSTD,skip=0):
    TranslateGold = MatchToGoldStd(FileCompare,FileGoldSTD)
    GoldSTD = pd.read_csv(FileGoldSTD)
    FateTrack = pd.read_csv(FileCompare)

    FateTrack = FateTrack.loc[FateTrack["frameNumber"]==frame,]
    FateTrack = FateTrack.reset_index(drop=True)
    GoldSTD = GoldSTD.loc[GoldSTD["frameNumber"]==frame,]
    GoldSTD = GoldSTD.reset_index(drop=True)

    correct=0
    incorrect=0
    for obj in range(FateTrack.shape[0]):
        FateID = FateTrack.loc[obj,"pointID"]
        FateParent = FateTrack.loc[obj,"parentID"]
        transGoldID = TranslateGold.loc[TranslateGold["FateID"]==FateID,"GoldID"].values[0] ;
        transGoldParent = TranslateGold.loc[TranslateGold["FateID"]==FateParent,"GoldID"] ;
        if not(transGoldParent.empty):
            transGoldParent = transGoldParent.values[0]
            actualGoldParent = GoldSTD.loc[GoldSTD["pointID"] == transGoldID,"parentID"]
            if (not(actualGoldParent.empty | math.isnan(actualGoldParent.values[0]))):
                actualGoldParent = int(actualGoldParent.values[0])
                if(actualGoldParent == transGoldParent):
                    correct = correct+1
                else:
                    incorrect = incorrect+1
    results = pd.DataFrame([frame, skip, correct, incorrect]).T
    results = results.rename(columns={0:"Frame",1:"Skip",2:"Correct",3:"Incorrect"})
    return(results)

def AssembleAccMeasurements(FileCompare,FileGoldSTD,skip=0):
    GoldSTD = pd.read_csv(FileGoldSTD)
    maxFrame = np.max(GoldSTD["frameNumber"])

    completeResults = pd.DataFrame()
    for frame in (np.array(range(1,maxFrame))+1):
        tmpFrame = CheckAccuracy(frame=frame,FileCompare=FileCompare,FileGoldSTD=FileGoldSTD,skip=skip)
        completeResults = completeResults.append(tmpFrame)
    completeResults = completeResults.reset_index(drop=True)
    return(completeResults)

def redefineGold(FileGoldSTD, outfilename, skip = 1,startTime = 0):
    GoldSTD = pd.read_csv(FileGoldSTD)

    sub = startTime+1
    maxFrame = np.max(GoldSTD['frameNumber'])
    frames_to_keep = np.array(range(startTime+1,maxFrame+1,skip+1))

    starter_frame = frames_to_keep[0]
    other_frames = frames_to_keep[1:]

    newGoldSTD = GoldSTD.loc[GoldSTD["frameNumber"].isin(other_frames),:]
    newGoldSTD = newGoldSTD.reset_index(drop=True)
    starterGold = GoldSTD.loc[GoldSTD["frameNumber"]==starter_frame,:]
    starterGold = starterGold.reset_index(drop=True)
    starterGold["parentID"] = "NaN"

    pointsNew = pd.concat([starterGold, newGoldSTD])["pointID"].values
    framesOld =  np.unique(newGoldSTD["frameNumber"])

    transmitFrame = pd.DataFrame()

    for i in range(newGoldSTD.shape[0]):
        tmpID = newGoldSTD.loc[i,"pointID"]
        tmpParent = GoldSTD.loc[GoldSTD["pointID"] ==tmpID,"parentID"].values[0]
        tmpFrame = newGoldSTD.loc[i,"frameNumber"]
        tmpX = newGoldSTD.loc[i,"xCoord"]
        tmpY = newGoldSTD.loc[i,"yCoord"]

        while (not((tmpParent in pointsNew)|(np.isnan(tmpParent)))):
            tmpParent = GoldSTD.loc[GoldSTD["pointID"] ==tmpParent,"parentID"].values[0]

        transmitFrame = transmitFrame.append(pd.DataFrame([tmpID,tmpFrame,tmpX,tmpY,tmpParent]).T)

    transmitFrame = transmitFrame.rename(columns={0:"pointID",1:"frameNumber",2:"xCoord",3:"yCoord",4:"parentID"})
    transmitFrame["annotation"] = "none"
    transmitFrame = transmitFrame.reset_index(drop=True)

    newFrame = pd.concat([starterGold,transmitFrame])
    newFrame = newFrame.reset_index(drop=True)
    newFrame['oldFrame'] = newFrame['frameNumber']
    newFrame = newFrame.drop(['frameNumber'],axis=1)

    changeFrame = pd.DataFrame(data={'oldFrame': frames_to_keep, 'frameNumber': frames_to_keep.argsort()+1})
    FinalFrame = pd.merge(newFrame,changeFrame)
    FinalFrame = FinalFrame.loc[:,["pointID","frameNumber","xCoord","yCoord","parentID","annotation"]]
    FinalFrame["pointID"] = FinalFrame["pointID"].astype('int64')

    FinalFrame.to_csv(outfilename,index=False)

def getParent(time,mainpath="./"):

    # Generate past features
    frame_start_file = mainpath+str(time)+'_seg.npy'
    frame_end_file = mainpath+str(time+1)+'_seg.npy'
    feature_start_file = mainpath+str(time)+'.feature.csv'
    feature_end_file = mainpath+str(time+1)+'.feature.csv'



    propies = generateLinks(frame_start_file, frame_end_file,timepoint=time,nnDist = 4,costMax=35,pathtoimage=mainpath)
    tmpdivs = DivisionCanditates(propMtx=propies,filename_t0=frame_start_file,MDAR_thresh=0.75, SDis_thresh=30,mS_Disp=1,
                                    filename_t1=frame_end_file,timepoint=time,pathtoimage=mainpath)
    finaldivs = UpdateConnectionsDiv(propies, tmpdivs)

    translation_table = TranslationTable(frame_start_file, frame_end_file,DivisionTable=finaldivs, timepoint=time,pathtoimage=mainpath)
    minCost_table = SolveMinCostTable(frame_start_file, frame_end_file,DivisionTable=finaldivs, timepoint=time,pathtoimage=mainpath)
    finTable = ReviewCostTable(minCostFlowtable = minCost_table, timepoint=time)
    masterConnects = TranslateConnections(finTable, translation_table, timepoint=time, preference="Master_ID")

    return(masterConnects)
