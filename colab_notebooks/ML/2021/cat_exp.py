# concatenate feature files from different frames, same experiment, into one csv
# concatenated file is in data folder, original files are in Trident folder

import os, sys
import glob
import pandas as pd

data_dir = sys.argv[1] # ~/Desktop/ml_data/data_0527_9dil/
output_dir = sys.argv[2] # ~/Desktop/ml_data/data_0527_9dil/
sample = sys.argv[3] # B3, B4
n_samples = int(sys.argv[4]) # 5, 7


os.chdir(data_dir)

filetypes = ['nuclei.feature', 'yfp_measurements', 'cy5_measurements', 'cy3_measurements']
for j in filetypes:
    all_filenames = []
    for k in range(1, n_samples):
        k_str = str(k)
        all_filenames.append(data_dir + sample + '_' + k_str + '/' + sample + '_' + k_str + '_' + j + '.csv')
    all_filenames_df = pd.DataFrame(all_filenames)
    all_filenames_df.to_csv('~/Desktop/output.csv')
    #combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
    #export to csv
    combined_csv.to_csv(output_dir + sample + '_' + j + ".csv", index=False)

##filetypes = ['nuclei', 'gfp', 'cy5']
##for j in filetypes:
##    all_filenames = []
##    for k in range(1, n_samples):
##        k_str = str(k)
##        if j == 'nuclei':
##            all_filenames.append(data_dir + j + '_' + sample + '_00' + k_str + '.feature.csv')
##        else:
##            all_filenames.append(data_dir + j + '_' + sample + '_00' + k_str + '.intensity.csv')
##    all_filenames_df = pd.DataFrame(all_filenames)
##    all_filenames_df.to_csv('~/Desktop/output.csv')
##    #combine all files in the list
##    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])
##    #export to csv
##    combined_csv.to_csv(output_dir + sample + '_' + j + ".csv", index=False)

