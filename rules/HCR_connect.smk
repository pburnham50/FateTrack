### HCR_connect.smk
from bin.HCR_connect import *

# Connect frames

rule HCR_connect:
  input:
#    TLmask='{path}/segment/{sample}/'+str(config['correction']['endFrame']-1)+'_seg.npy',
    TLcheck='{path}/segment/{sample}/0.tif',
    HCRmask='{path}/HCR/{sample}_finalState_Nuclear_seg.npy'
  output:
    HCRsub='{path}/HCR/{sample}_HCR_connections_MasterID.csv'
  params:
    timepoint=(config['correction']['endFrame']-1),
    nnDist = config['connection']['nnDist'],
    costMax = config['connection']['costMax'],
    mN_Int = config['connection']['mN_Intensity'],
    mN_Ecc = config['connection']['mN_Eccentricity'],
    mN_Area = config['connection']['mN_Area'],
    mN_Disp = config['connection']['mN_Displacement'],
    mS_Area = config['connection']['mS_Area'],
    mS_Ecc = config['connection']['mS_Eccentricity'],
    mS_Int = config['connection']['mS_Intensity'],
    mS_Disp = config['connection']['mS_Displacement'],
    MDAR_thresh = config['connection']['MDAR_thresh'],
    SDis_thresh = config['connection']['SDis_thresh'],
    openingCost = config['connection']['openingCost'],
    closingCost = config['connection']['closingCost']
  run:
    HCR_connect(sampleName=wildcards.sample, TLlast_mask=wildcards.path+'/segment/'+wildcards.sample+'/'+str(params.timepoint)+'_seg.npy',
                HCR_mask=input.HCRmask, timepoint=params.timepoint,
                nnDist = params.nnDist, costMax = params.costMax, mN_Int = params.mN_Int, mN_Ecc = params.mN_Ecc,
                mN_Area = params.mN_Area, mN_Disp = params.mN_Disp, mS_Area = params.mS_Area, mS_Ecc = params.mS_Ecc,
                mS_Int = params.mS_Int, mS_Disp = params.mS_Disp, MDAR_thresh = params.MDAR_thresh,
                SDis_thresh = params.SDis_thresh, openingCost = params.openingCost, closingCost = params.closingCost)
