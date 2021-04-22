### segment.smk

# Breaks single TIFF stack into individual files

rule segment:
  input:
    get_frame_files("tmp/GFP20mtest.frame.info.csv", 'frame_image')
  output:
    get_frame_files("tmp/GFP20mtest.frame.info.csv", 'frame_seg')
  params:
    size=19.0,
    name=config['exp']
  shell:
    "python -m cellpose --dir results/{params.name}/frames/ --pretrained_model nuclei --diameter {params.size} ;"
