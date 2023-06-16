# concatenate feature files from different frames, same experiment, into one csv
# concatenated file is in data folder, original files are in Trident folder

import os, sys
import glob
import pandas as pd

data_dir = sys.argv[1] # ~/Desktop/ml_data/HCR_0610/tables/
output_dir = sys.argv[2] # ~/Desktop/ml_data/HCR_0610/507_D3_1_w1/
sample = sys.argv[3] # sub
n_samples = int(sys.argv[4]) # 16


os.chdir(data_dir)

filetypes = ['1_nuclei.feature', '2_nuclei.feature', '1YFP_measurements', '1CY5_measurements', '1CY3_measurements', '2YFP_measurements', '2CY5_measurements', '2CY3_measurements']
for j in filetypes:
    combined_csv = pd.DataFrame()
    for k in range(0, n_samples):
        k_str = str(k)
        filename = data_dir + sample + '_' + k_str + '_' + j + '.csv'
        f_csv = pd.read_csv(filename)
        l = len(f_csv.columns)
        #subimage_col = [k] * l
        #f_csv['sample'] = sample_col
        f_csv.insert(0, "subimage", k)
        combined_csv = pd.concat([combined_csv, f_csv])
    #export to csv
    combined_csv.to_csv(output_dir + j + ".csv", index=False)

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

