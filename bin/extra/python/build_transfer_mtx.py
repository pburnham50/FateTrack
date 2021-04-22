# prints csv file which is used

import os
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser('Generate file of image properties')

    parser.add_argument('--dir', type = str, dest='direct', help = 'Directory containing tiffs.')
    parser.add_argument('--name', type = str, dest='out',help = 'Name to output.')

    args = parser.parse_args()


    numImages = int(os.popen('ls '+args.direct+'*.tif | wc -l').read())

    connectFrame = pd.DataFrame()
    for i in range(1,numImages):
        tmpFrame=pd.DataFrame([[str(i-1),i-1,i,args.direct+str(i-1)+".tif",args.direct+str(i)+".tif",
                                args.direct+str(i-1)+"_seg.npy",args.direct+str(i)+"_seg.npy",
                                args.direct+str(i-1)+".knex.csv"]],
                                columns=('id','frame_start','frame_end','frame_image_start',
                                        'frame_image_end','frame_seg_start','frame_seg_end',
                                        'connection_file'))
        connectFrame = connectFrame.append(tmpFrame)

    #outframe = connectFrame[['frame']]
    connectFrame.to_csv("tmp/"+args.out+".connect.info.csv", index = False, header=True)

    onlyFrame = pd.DataFrame()
    for i in range(0,numImages):
        tmpFrame=pd.DataFrame([[i,args.direct+str(i)+".tif",args.direct+str(i)+"_seg.npy"]],
                                columns=('frame','frame_image','frame_seg'))
        onlyFrame = onlyFrame.append(tmpFrame)

    onlyFrame.to_csv("tmp/"+args.out+".frame.info.csv", index = False, header=True)

    return

if __name__ == '__main__':
    main()
