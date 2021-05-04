#FT_connect
from FT_connect_functions import *
from FT_connect_config import *

for i in range(end_frame):
    t0,t1,ffC = fullPipeline(file_1=str(i)+"_seg.npy",
                       timept_1=str(i),
                       file_2=str(int(i)+1)+"_seg.npy",
                       timept_2=str(int(i)+1))

    rlg = RLGwrap(finalConnections=ffC,FeatureFrame_t0=t0,FeatureFrame_t1=t1,time0=int(i))

    if not os.path.exists(rlg_out_dir):
        os.makedirs(rlg_out_dir)

    rlg.to_csv(rlg_out_dir+"out_"+str(i)+".csv",index=False)

appended_data = pd.DataFrame()
frameTable = pd.DataFrame()
for i in range(end_frame):
    data = pd.read_csv(rlg_out_dir+"out_"+str(i)+".csv")
    appended_data = appended_data.append(data)
    frameTable = frameTable.append(pd.DataFrame(['cy5_time'+str(i)+'.tif',str(i),'cy5',str(i+1)]).T)
for i in range(end_frame+1):
    copyfile(pathToSegs+str(i)+'.tif',rlg_out_dir+'cy5_time'+str(i)+'.tif')


frameTable = frameTable.rename(columns={0:'fileName',1:'time',2:'wavelength',3:'frameNumber'})
frameTable.to_csv(rlg_out_dir+"fileTable.csv",index=False)
appended_data.to_csv(rlg_out_dir+"out.csv",index=False)
