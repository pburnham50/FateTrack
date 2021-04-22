### breakup.smk
from bin.fatetrack_segment import *

# Breaks single TIFF stack into individual files

rule TL_segment:
  input:
    nuclear='{path}/correction/{sample}_Nuclear_corr.tif'
  output:
    nuclearCorr='{path}/segment/{sample}/0.tif',
    segmentDirectory=directory('{path}/segment/{sample}/')
  params:
    flowThreshold=config['segment']['flowThreshold'],
    cellProbabilityThreshold=config['segment']['cellProbabilityThreshold'],
    nucDiam=config['segment']['nuclearDiameter']
  run:
    outpath = str(wildcards.path)+'/segment/'+str(wildcards.sample)
    segmentCellpose(inputImage=input.nuclear, cellpose_out=outpath, nucDiameter=params.nucDiam, cellprob_threshold=params.cellProbabilityThreshold, flow_threshold=params.flowThreshold)
