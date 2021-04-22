### HCR_subImage.smk
from bin.HCR_matching import *

# Connect frames

rule HCR_subImage:
  input:
    TLimage='{path}/correction/{sample}_Nuclear_corr.tif',
    HCRimage=config['experiment']['HCR_A594']
  output:
    HCRsub='{path}/HCR/{sample}_finalState_Nuclear.tif',
    HCRcoords='{path}/HCR/{sample}_HCRsubCoords.txt',
    coarseHCRcoords='{path}/HCR/{sample}_HCRsubCoords_coarse.txt',
    redistrictHCRcoords='{path}/HCR/{sample}_HCRsubCoords_redistrict.txt',
    HCRmask='{path}/HCR/{sample}_finalState_Nuclear_seg.npy'
  params:
    endFrame=config['correction']['endFrame'],
    stitch_tiles=config['HCR']['stitchTiles'],
    penalty_multiplier=config['HCR']['penaltyCoeff'],
    cost_inclusion = config['HCR']['costInclusion'],
    flowThreshold=config['segment']['flowThreshold'],
    cellProbabilityThreshold=config['segment']['cellProbabilityThreshold'],
    rows=config['experiment']['TLrows'],
    cols=config['experiment']['TLcolumns'],
    pxO=config['experiment']['pxOverlap'],
    nucDiam=config['segment']['nuclearDiameter']
  run:
    getHCRsubImage(outFile=output.HCRsub, TLnuclearImage=input.TLimage, \
                    Coarse_HCRsubcoordsFile=output.coarseHCRcoords, Restrict_HCRsubcoordsFile=output.redistrictHCRcoords,
                    HCRnuclearImage=input.HCRimage,frame = (params.endFrame-1), \
                    Fine_HCRsubcoordsFile=output.HCRcoords, stitch_tiles=params.stitch_tiles, \
                    penalty_multiplier=params.penalty_multiplier, \
                    cost_inclusion = params.cost_inclusion, totalRows=params.rows, \
                    totalCols=params.cols, pixelOverlap=params.pxO) ;

    getHCRmask(HCRsubImage=output.HCRsub, nucDiameter=params.nucDiam, \
                cellprob_threshold=params.cellProbabilityThreshold, \
                flow_threshold=params.flowThreshold) ;
