# FateTrack
Tools to track and classify cell states in organoid monolayers.

![GitHub Logo](/design/FateTrack_logo-11.png)

## Summary of pipeline
1. Stabilization and Denoising
2. Segmentation 
3. Feature Extraction
4. Connection
5. Classification - Static and dynamic features are compared to training data to determine cell-type.

## Requirements
* TIFF stack of organoid monolayer in nuclear fluorescent channel (stack through time) and birghtfield or phase-contrast.
* metadata file including time infomration

## Procedure

1. Create conda environment:
```
$ conda env create -f environment.yml
```

2. Edit configuration file to include relevant information:
