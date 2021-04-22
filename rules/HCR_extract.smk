### HCR_extract.smk
from bin.HCR_extract import *
import pandas as pd

rule HCR_extract:
  input:
    HCRcoords='{path}/HCR/{sample}_HCRsubCoords.txt',
    HCRmask='{path}/HCR/{sample}_finalState_Nuclear_seg.npy',
    HCRimage=config['experiment']['HCR_A594'],
    connectionPath=directory('{path}/connect/{sample}/'),
    finalConnect='{path}/HCR/{sample}_HCR_connections_MasterID.csv'
  output:
    HCRmeasures='{path}/HCR/{sample}_finalState_meaurements.csv',
    VoroMask='{path}/HCR/{sample}_Voro_seg.npy',
    VoroImage='{path}/HCR/{sample}_Voro.png',
    VoroTIF='{path}/HCR/{sample}_Voro_final.tif',
    VoroTrans='{path}/HCR/{sample}_Clps2Voro.csv',
    allConnect='{path}/HCR/{sample}_totalConnections.csv'
  params:
    timepoint=(config['correction']['endFrame']),
    channelList = ['DAPI','YFP','A594','CY3', 'CY5'],
    dilationList = [1], # for expanding masks to measure intensity.
    MaxVoroArea = 1600
  run:
    intDF = assembleFinalState(HCR_mask_file=input.HCRmask,
                                HCR_image_file=input.HCRimage,
                                HCR_coords_file=input.HCRcoords,
                                Voro_mask_file=output.VoroMask,
                                Voro_image_file=output.VoroImage,
                                Voro_image_file_final=output.VoroTIF,
                                Voro_Transfer_file=output.VoroTrans,
                                channelList=params.channelList,
                                final_timepoint=params.timepoint,
                                dilations=params.dilationList,
                                MaxVoroArea=params.MaxVoroArea)
    print(input.connectionPath, input.finalConnect, params.timepoint,wildcards.sample)
    connex = getAllConnections(pathToConnects=input.connectionPath+'/', #'results/'++'/connect/'+wildcards.sample+'/',
                                finalConnection=input.finalConnect,
                                final_timepoint=params.timepoint,
                                sampleName=wildcards.sample)

    intDF.to_csv(output.HCRmeasures)
    connex.to_csv(output.allConnect)
