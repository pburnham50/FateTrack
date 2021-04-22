### breakup.smk

# Breaks single TIFF stack into individual files
storage = config['storage'] #'/Volumes/PSB_2020_2/FateTrack/select_data/timelapse/FT_317_Day4_1_w2'

rule TL_correct:
  input:
    nuclear=storage+'/{sample}_A594.tif',
    bf=storage+'/{sample}_Brightfield.tif'
  output:
    nuclearCorr='{path}/correction/{sample}_Nuclear_corr.tif',
    bfCorr='{path}/correction/{sample}_Brightfield_corr.tif'
  params:
    bkgrd=config['correction']['backgroundImage'],
    stFrame=config['correction']['startFrame'],
    endFrame=config['correction']['endFrame'],
    pxCrop=config['correction']['pixelCroppage']
  shell:
    " python bin/FT_correct.py --nucImage {input.nuclear} --bfImage {input.bf} --outpath results/{wildcards.sample}/correction/ --backgroundImage {params.bkgrd} --cropPlus {params.pxCrop} --startFrame {params.stFrame} --stopFrame {params.endFrame} ; "
