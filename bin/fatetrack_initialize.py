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
import fatetrack_connections
import fatetrack_review
from fatetrack_connections import buildFeatureFrame, buildOffsetFrame, generateCandidates, generateLinks, DivSimScore, DivSetupScore, DivisionCanditates, UpdateConnectionsDiv, TranslationTable, SolveMinCostTable, ReviewCostTable
from fatetrack_review import TranslateConnections, RajTLG_wrap, MatchToGoldStd, CheckAccuracy, AssembleAccMeasurements, redefineGold, getParent
from shutil import copyfile

################################################################################################################################
################ MAIN
################################################################################################################################

def main():
    parser = argparse.ArgumentParser('Generate a connection table that can then be manually reviewed using Time Lapse GUI.')

    parser.add_argument('--path', type = str, dest='mainpath', help = 'Path to segmentation frames.')
    parser.add_argument('--outpath', type = str, dest='outpath', help = 'Path to connection info.')
    parser.add_argument('--maxTime', type = int, dest='maxtime',help = 'The latest frame number, as an integer (e.g. \'10\')')
    parser.add_argument('--nearestNeighborDistance', type = int, default=35, dest='nearestNeighborDistance',help = 'Initalizes radius from each nucleus/cell to look. e.g. 4 = 4x(distance to nearest neighbor)')
    parser.add_argument('--linkageCostMaximum', type = int, default=35, dest='linkageCostMaximum',help = 'Sets a limit for links to keep by looking at energy.')
    parser.add_argument('--openingCost', type = int, default=30, dest='openCost', help = 'Cost to open connections without connection to a previous segmentation.')
    parser.add_argument('--closingCost', type = int, default=30, dest='closeCost', help = 'Cost to close connections without connection to a next segmentation.')
    parser.add_argument('--parentChildAreaChange', type = float, default=0.75, dest='parentChildAreaChange',help = 'The ratio of mother to daughter cells area, sets a limit for daughter-pairs to evaluate as division DivCandidatesMtx.')
    parser.add_argument('--siblingThreshold', type = int, default=30,dest='siblingThreshold',help = 'Maximum allowable sibling cost for two nuclei.')
    parser.add_argument('--multiplierDivisionDisplacement', type = int, default=1, dest='multiplierDivisionDisplacement',help = 'Weighs the differentiation score between potential siblings by their displacement.')
    parser.add_argument('--multiplierDivisionArea', type = int, default=10, dest='multiplierDivisionArea',help = 'Weighs the differentiation score between potential siblings by their area difference.')
    parser.add_argument('--multiplierDivisionEccentricity', type = int, default=2, dest='multiplierDivisionEccentricity',help = 'Weighs  the differentiation score between potential siblings by their eccentricity difference.')
    parser.add_argument('--multiplierDivisionIntensity', type = int, default=2, dest='multiplierDivisionIntensity',help = 'Weighs the differentiation score between potential siblings by their intensity difference.')
    parser.add_argument('--multiplierConnectionDisplacement', type = int, default=1, dest='multiplierConnectionDisplacement',help = 'Weighs the differentiation score between potential connections by their displacement.')
    parser.add_argument('--multiplierConnectionArea', type = int, default=25, dest='multiplierConnectionArea',help = 'Weighs the differentiation score between potential connections by their area difference.')
    parser.add_argument('--multiplierConnectionEccentricity', type = int, default=4, dest='multiplierConnectionEccentricity',help = 'Weighs the differentiation score between potential connections by their eccentricity difference.')
    parser.add_argument('--multiplierConnectionIntensity', type = int, default=10, dest='multiplierConnectionIntensity',help = 'Weighs the differentiation score between potential connections by their intensity difference.')

    args = parser.parse_args()

    mainpath = args.mainpath+'/'
    outpath = args.outpath+'/'
    maxtime = args.maxtime
    nearestNeighborDistance = args.nearestNeighborDistance
    linkageCostMaximum = args.linkageCostMaximum
    parentChildAreaChange = args.parentChildAreaChange
    siblingThreshold = args.siblingThreshold
    multiplierDivisionDisplacement = args.multiplierDivisionDisplacement
    multiplierDivisionArea = args.multiplierDivisionArea
    multiplierDivisionEccentricity = args.multiplierDivisionEccentricity
    multiplierDivisionIntensity = args.multiplierDivisionIntensity
    multiplierConnectionDisplacement = args.multiplierConnectionDisplacement
    multiplierConnectionArea = args.multiplierConnectionArea
    multiplierConnectionEccentricity = args.multiplierConnectionEccentricity
    multiplierConnectionIntensity = args.multiplierConnectionIntensity
    openCost = args.openCost
    closeCost = args.closeCost

    allTranslate = pd.DataFrame()
    allConX = pd.DataFrame()
    frameTable = pd.DataFrame()
    masterConnects_Master = pd.DataFrame()

    if not os.path.exists(outpath+'review/'):
        os.mkdir(outpath+'review/')

    for time in range(0,(maxtime+1)):

        # Generate past features
        frame_start_file = mainpath+str(time)+'_seg.npy'
        frame_end_file = mainpath+str(time+1)+'_seg.npy'

        propies = generateLinks(frame_start_file, frame_end_file,timepoint=time,nnDist = nearestNeighborDistance,costMax=linkageCostMaximum,pathtoimage=mainpath,mN_Int = multiplierConnectionIntensity, mN_Ecc=multiplierConnectionEccentricity, mN_Area=multiplierConnectionArea, mN_Disp=multiplierConnectionDisplacement)
        tmpdivs = DivisionCanditates(propMtx=propies,filename_t0=frame_start_file,MDAR_thresh=parentChildAreaChange, SDis_thresh=siblingThreshold,mS_Disp=multiplierDivisionDisplacement,mS_Area = multiplierDivisionArea, mS_Ecc = multiplierDivisionEccentricity, mS_Int = multiplierDivisionIntensity,filename_t1=frame_end_file,timepoint=time,pathtoimage=mainpath)
        finaldivs = UpdateConnectionsDiv(propies, tmpdivs)

        translation_table = TranslationTable(frame_start_file, frame_end_file,DivisionTable=finaldivs, timepoint=time,pathtoimage=mainpath)
        minCost_table = SolveMinCostTable(frame_start_file, frame_end_file,DivisionTable=finaldivs, timepoint=time,pathtoimage=mainpath,OpeningCost = openCost, ClosingCost = closeCost)
        finTable = ReviewCostTable(minCostFlowtable = minCost_table, timepoint=time)
        masterConnects = TranslateConnections(finTable, translation_table, timepoint=time, preference="RajTLG_ID")
        masterConnects_Master_tmp = TranslateConnections(finTable, translation_table, timepoint=time, preference="Master_ID")
        masterConnects_Master_tmp.to_csv(outpath+'Connections_MasterID_'+str(time)+'.csv', index=False)

        tmpRLG = RajTLG_wrap(filename_t0=frame_start_file, filename_t1=frame_end_file, timepoint=time,path=mainpath,ConnectionTable=finTable,TranslationTable=translation_table)
        allConX = allConX.append(tmpRLG)
        allTranslate = allTranslate.append(translation_table)
        copyfile(mainpath+str(time)+'.tif',outpath+'review/cy5_time'+str(time)+'.tif')
        tmpFrameInfo = pd.DataFrame(['cy5_time'+str(time)+'.tif',str(time),'cy5',str(time+1)]).T
        frameTable = frameTable.append(tmpFrameInfo)

    copyfile(mainpath+str(maxtime+1)+'.tif',outpath+'review/cy5_time'+str(time+1)+'.tif')
    frameTable = frameTable.append(pd.DataFrame(['cy5_time'+str(maxtime+1)+'.tif',str(maxtime+1),'cy5',str(maxtime+2)]).T)
    frameTable = frameTable.rename(columns={0:'fileName',1:'time',2:'wavelength',3:'frameNumber'})
    frameTable.to_csv(outpath+"review/fileTable.csv",index=False)
    allTranslate.to_csv(outpath+"translationTable.csv",index=False)
    allConX.to_csv(outpath+"review/out.csv",index=False)

    return

if __name__ == '__main__':
    main()
