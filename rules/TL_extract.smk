### TL_connect.smk
from bin.fatetrack_nucFeatureExtraction import *

# Connect frames

rule TL_extract:
  input:
    bfImage='{path}/correction/{sample}_Brightfield_corr.tif',
    segPath=directory('{path}/segment/{sample}/')
  output:
    features='{path}/features/{sample}.0_staticFeatures.csv'
  run:
    outpath = str(wildcards.path)+'/features/'
    staticFeatureExt(pathToSegmenations = input.segPath, imageFile = input.bfImage, outpath = outpath, sampleName = str(wildcards.sample))
