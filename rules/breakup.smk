### breakup.smk

# Breaks single TIFF stack into individual files

rule breakup_stack:
  input:
    "volumes/{expfov}.tif"
  output:
    singletif="results/{expfov}/frames/0.tif",
    mapping="tmp/frame_mapping.csv"
  shell:
    "magick convert {input} results/{wildcards.expfov}/frames/%d.tif ; "
    "python bin/extra/python/build_transfer_mtx.py --dir results/{wildcards.expfov}/frames/ --csv {output.mapping} ;"
