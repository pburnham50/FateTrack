### TL_connect.smk

# Connect frames

rule TL_connect:
  input:
    segmentDirectory=directory('{path}/segment/{sample}/'),
    firstImage='{path}/segment/{sample}/0.tif'
  output:
    connectDirectory=directory('{path}/connect/{sample}/'),
    reviewDirectory=directory('{path}/connect/{sample}/review/'),
    translationTable='{path}/connect/{sample}/translationTable.csv'
  params:
    frameCount=(config['correction']['endFrame']-2)
  shell:
    " python bin/fatetrack_initialize.py --path {input.segmentDirectory} --outpath {output.connectDirectory} --maxTime {params.frameCount} ; "
