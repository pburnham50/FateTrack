Here is a copy of all Colab notebooks and supplemental models used for the FateTrack pipeline.

To run the pipeline, use
FateTrack_Processing_Current.ipynb
This script has dependencies of (each is contained within its own folder):
Trained CARE model: in CARE_model
Trained Cellpose model: in cellpose_model
Trained FateTrack connection algorithm (min-cost-flow): in connection_scripts

Additional scripts were created to create files for annotating data. These are located in the annotation folder and are:
annotation.ipynb: generated annotation files in 2022
annotation_2023.ipynb: generated annotation files in 2023
createRLGfile.ipynb: create RajLabGUI file format from FateTrackObject to use with timelapseGUI
FateTrack_Processing_Manual_Annotations.ipynb: adaptation of FateTrack_Processing_Current.ipynb that reads in the from manual annotation files (out.csv created with timelapseGUI or FateTrackSelector) and incorporates into the FateTrackObject

retired_processing_scripts: old versions of FateTrack_Processing_Current.ipynb
ML: any notebooks and scripts containing ML models or to process outputs from ML models
AnnData: adaptation of FateTrack_Processing_Current.ipynb into the AnnData file format
