library(ggtern)
library(ggplot2)
library(umap)

setwd("/Users/pburnham/Documents/GitHub/FateTracking/")

fileList = list.files(path="analysis/finalStateHCR_v4/",pattern = "*.csv")
blacklistFiles = c(97,116,117,157:169)


checking = c()
for (i in fileList){
  tmp = data.frame(read.csv(paste0("analysis/finalStateHCR_v4/",i)))
  tmp$sample = paste(strsplit(i,"_")[[1]][1:6],collapse = '_')
  tmp$cellID = paste0(tmp$sample,"|",tmp$Master_ID)
  colnames(tmp)[35] = "Master_ID_End"
  checking = rbind(checking,tmp)
}

rangeChange <- function(x){(x-min(x))/(max(x)-min(x))}

checking$CY3 = (checking$CY3_sum_intensity_1)#/min(checking$CY3_background))
checking$CY5 = (checking$CY5_sum_intensity_1)#/min(checking$CY5_background))
checking$YFP = (checking$YFP_sum_intensity_1)#/min(checking$YFP_background))
checking$CY3 = (checking$CY3/min(checking$CY3+1))
checking$CY5 = (checking$CY5/min(checking$CY5+1))
checking$YFP = (checking$YFP/min(checking$YFP+1))
checking$CY3norm = checking$CY3 / (checking$CY3 + checking$CY5 + checking$YFP)
checking$CY5norm = checking$CY5 / (checking$CY3 + checking$CY5 + checking$YFP)
checking$YFPnorm = checking$YFP / (checking$CY3 + checking$CY5 + checking$YFP)

checking$magnitude =  sqrt((checking$YFP)**2 + (checking$CY5)**2 + (checking$CY3)**2)
checking_high = checking[checking$magnitude>quantile(checking$magnitude,probs = .75),]
checking_low = checking[checking$magnitude<quantile(checking$magnitude,probs = .05),]

checking_high = checking_high[,(ncol(checking_high)-9):(ncol(checking_high))]
checking_low = checking_low[,(ncol(checking_low)-9):(ncol(checking_low))]
checking = checking[,(ncol(checking)-9):(ncol(checking))]

#ggtern(checking_high,aes(CY3,CY5,YFP)) +
#  geom_point(size=.5) + 
#  stat_density_tern(geom='polygon',
#                    aes(col=..level..),fill='transparent') +
#  scale_color_gradient(low='grey',high='red') +
#  limit_tern(1,1,1)+
#  theme_dark() +
#  theme_showarrows() +
#  theme_clockwise()+
#  theme(legend.position = "none")


goblet = (checking_high[(checking_high$CY5norm > .5)&(checking_high$YFPnorm <.3)&(checking_high$CY3norm <.3),])
entero = (checking_high[(checking_high$YFPnorm > .8)&(checking_high$CY5norm <.1)&(checking_high$CY3norm <.1),])
paneth = (checking_high[(checking_high$CY3norm > .5)&(checking_high$CY5norm <.3)&(checking_high$YFPnorm <.3),])
GEmix = (checking_high[(checking_high$YFPnorm > .35)&(checking_high$YFPnorm < .6)&
                        (checking_high$CY5norm >.35)&(checking_high$CY5norm <.6),])
GPmix = (checking_high[(checking_high$CY3norm > .35)&(checking_high$CY3norm < .6)&
                         (checking_high$CY5norm >.35)&(checking_high$CY5norm <.6),])

minAmount = 2*min(c(nrow(goblet),nrow(entero),nrow(GEmix)))



ggplot(checking,aes(CY5,YFP))+geom_point(size=0.8,alpha=.5)+
  theme_bw()
ggplot(checking,aes(CY5norm,YFPnorm))+
  geom_point(size=0.8,alpha=.5,aes(size=CY3norm))+
  geom_point(data=entero,size=0.8,alpha=.5,col="green",aes(size=CY3norm))+
  geom_point(data=goblet,size=0.8,alpha=.5,col="red",aes(size=CY3norm))+
  geom_point(data=paneth,size=0.8,alpha=.5,col="blue",aes(size=CY3norm))+
  geom_point(data=checking_low,size=0.8,alpha=.5,col="gold",aes(size=CY3norm))+
  theme_dark()
######

getStaticFeatures = function(Master_ID,experiment,track,finalID){
  cell_tmp_id = Master_ID
  frame = strsplit(cell_tmp_id,"_")[[1]][1]
  obj = strsplit(cell_tmp_id,"_")[[1]][2]
  
  cellFeaturesDir = paste("results",experiment,"features","",sep = "/")
  if(file.exists(paste0(cellFeaturesDir,experiment,".",frame,"_staticFeatures.csv"))){
    staticFeatures = data.frame(read.csv(paste0(cellFeaturesDir,experiment,".",frame,"_staticFeatures.csv")))
    tmp_feature = staticFeatures[staticFeatures$label == obj,]
    tmp_feature$trackID = track # as.character(TL_HCR_connect$TrackID)
    tmp_feature$finalID = finalID
    tmp_feature$sample = experiment
    return(tmp_feature)
  }
}


#####

gobs = c()
for (i in 1:nrow(goblet)){
  sampleGob = goblet[i,]
  CellID  = as.character(as.matrix(sampleGob$Master_ID_13))
  experiment =  as.character(as.matrix(sampleGob$sample))
  
  HCRDir = paste("results",experiment,"HCR","",sep = "/")
  TL_HCR_connect = data.frame(read.csv(paste0(HCRDir,experiment,"_totalConnections.csv")))
  TL_HCR_connect = TL_HCR_connect[TL_HCR_connect$Master_ID_13 == CellID,]
  if (nrow(TL_HCR_connect)>0){
    TL_characteristics  = as.matrix(TL_HCR_connect)[2:14]
  
    for(ii in 1:13){
      gobs = rbind(gobs, getStaticFeatures(TL_characteristics[ii],experiment = experiment, track =  as.character(TL_HCR_connect$TrackID), finalID=as.matrix(TL_HCR_connect)[ncol(TL_HCR_connect)]))
    }
  }
}


ents = c()
for (i in 1:nrow(entero)){
  sampleGob = entero[i,]
  CellID  = as.character(as.matrix(sampleGob$Master_ID_13))
  experiment =  as.character(as.matrix(sampleGob$sample))
  
  HCRDir = paste("results",experiment,"HCR","",sep = "/")
  TL_HCR_connect = data.frame(read.csv(paste0(HCRDir,experiment,"_totalConnections.csv")))
  TL_HCR_connect = TL_HCR_connect[TL_HCR_connect$Master_ID_13 == CellID,]
  if (nrow(TL_HCR_connect)>0){
    TL_characteristics  = as.matrix(TL_HCR_connect)[2:14]
    
    for(ii in 1:13){
      ents = rbind(ents, getStaticFeatures(TL_characteristics[ii],experiment = experiment, track =  as.character(TL_HCR_connect$TrackID), finalID=as.matrix(TL_HCR_connect)[ncol(TL_HCR_connect)]))
    }
  }
}



colnames(goblet)[1] = "finalID"
colnames(entero)[1] = "finalID"

goblet_2 = merge(goblet,gobs,by="finalID")
entero_2 = merge(entero,ents,by="finalID")

goblet_2$type = 'goblet'
entero_2$type = 'entero'
all = rbind(goblet_2,entero_2)
tracks = unique(all$trackID)
blacklistTracks = c()
for(j in 1:length(tracks)){
  tmptrack = (all[all$trackID == tracks[j],])
  timepoints = length(tmptrack$frame)
  if(any(sqrt((tmptrack$centroid.0[-1]-tmptrack$centroid.0[-timepoints])**2 + (tmptrack$centroid.1[-1]-tmptrack$centroid.1[-timepoints])**2 )>=15)){
    blacklistTracks = c(blacklistTracks,tracks[j])
  }
}
blacklistTracks = c(blacklistTracks,unique(all[(all$area<40)|(all$area>550),]$trackID))


remov.cols = c("frame" ,"centroid.0","centroid.1","trackID","type","intensity_image","label" , "bbox.0","bbox.1","bbox.2","bbox.3",'finalID',"sample","cellID","CY3","CY5",
               "YFP","CY3norm","CY5norm","YFPnorm","magnitude")

all = all[!(all$trackID %in% blacklistTracks),]
alltest = all[,!(colnames(all)%in%remov.cols)]
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]
alltest[is.na(alltest)]=0
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]

alltestscale = scale(alltest)
track_scales = cbind(alltestscale,all$trackID)
extrema = c()
for ( clmn in 1:(ncol(track_scales)-1)){
  list_extremes = which(abs(as.numeric(track_scales[,clmn]))>3.5)
  if(!is.null(list_extremes)){
    extrema = c(extrema, as.character(track_scales[list_extremes,ncol(track_scales)]))
  }
}
blacklist2 = unique(extrema)

all = all[!(all$trackID %in% blacklist2),]
alltest = all[,!(colnames(all)%in%remov.cols)]
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]
alltest[is.na(alltest)]=0
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]

alltestscale = scale(alltest)

#alltestscale = alltestscale[,- as.numeric(which(apply(alltestscale, 2, var) == 0))]


pca_alltest = prcomp((alltestscale))












##############
goblet$type = "goblet"
entero$type = "entero"
paneth$type = "paneth"

checking_low$type = "latent"
cellfro = rbind(entero[order(entero$magnitude,decreasing = T),][1:minAmount,],
                checking_low[order(checking_low$magnitude,decreasing = T),][1:minAmount,],
                goblet)
colnames(cellfro)[1] = "finalID"

allcells = c()
for (i in 1:nrow(cellfro)){
  sampleGob = cellfro[i,]
  CellID  = as.character(as.matrix(sampleGob$finalID))
  experiment =  as.character(as.matrix(sampleGob$sample))
  
  HCRDir = paste("results",experiment,"HCR","",sep = "/")
  TL_HCR_connect = data.frame(read.csv(paste0(HCRDir,experiment,"_totalConnections.csv")))
  TL_HCR_connect = TL_HCR_connect[ TL_HCR_connect[,ncol(TL_HCR_connect)] == CellID,]
  
  if (nrow(TL_HCR_connect)>0){
    print(paste0("good - ",experiment, " - ", sampleGob$type ))
    TL_characteristics  = as.matrix(TL_HCR_connect)[2:(ncol(TL_HCR_connect)-1)]
      
    for(ii in 1:13){
      allcells = rbind(allcells, getStaticFeatures(TL_characteristics[ii],experiment = experiment, track =  as.character(TL_HCR_connect$TrackID), finalID=as.matrix(TL_HCR_connect)[ncol(TL_HCR_connect)]))
    }
  }else{
    print(paste0("bad - ",experiment, " - ", sampleGob$type ))
  }
}


checking_2 = merge(cellfro,allcells,by=c("finalID","sample"))

check_export = checking_2#[,colnames(checking_2) %in% c('sample','finalID','CY3','CY5','YFP','CY3norm','CY5norm','YFPnorm','magnitutde','frame','label','trackID',"centroid.0","centroid.1","major_axis_length","type")]
check_export$MasterID = paste0(check_export$frame,"_",check_export$label)
names(check_export)[grep("centroid",colnames(check_export))] <- c('centroid0','centroid1')
check_export$fov = as.numeric(lapply(as.character(as.matrix(check_export$sample)),function(x) strsplit(x,split = 'Y')[[1]][2]))
check_export = check_export[!(check_export$fov %in% blacklistFiles),]


tracks = unique(check_export$trackID)
blacklistTracks = c()
for(j in 1:length(tracks)){
  tmptrack = (check_export[check_export$trackID == tracks[j],])
  timepoints = length(tmptrack$frame)
  if(any(sqrt((tmptrack$centroid0[-1]-tmptrack$centroid0[-timepoints])**2 + (tmptrack$centroid1[-1]-tmptrack$centroid1[-timepoints])**2 )>=(2*check_export$equivalent_diameter))){
    blacklistTracks = c(blacklistTracks,tracks[j])
  }
}


blacklistTracks = c(blacklistTracks,unique(check_export[(check_export$area<30)|(check_export$area>1000),]$trackID))


remov.cols = c("frame" ,"centroid0","centroid1","trackID","type","intensity_image","label" , "bbox.0","bbox.1","bbox.2","bbox.3",'finalID',"sample","cellID","CY3","CY5",
               "YFP","CY3norm","CY5norm","YFPnorm","magnitude","MasterID", "fov")
check_export =  check_export[!(check_export$trackID %in% blacklistTracks),]
check_export = check_export[,colnames(check_export) %in% c('sample','finalID','CY3','CY5','YFP','CY3norm','CY5norm','YFPnorm','magnitutde','frame','label','trackID',"centroid0","centroid1","major_axis_length","type")]

write.csv(x = check_export, "analysis/20210410_cellinfo_LGE.csv",quote = F,row.names = F,)

all = checking_2
#all = all[!(all$trackID %in% blacklistTracks),]
alltest = all[,!(colnames(all)%in%remov.cols)]
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]
alltest[is.na(alltest)]=0
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]

alltestscale = scale(alltest)
track_scales = cbind(alltestscale,all$trackID)
extrema = c()
for ( clmn in 1:(ncol(track_scales)-1)){
  list_extremes = which(abs(as.numeric(track_scales[,clmn]))>12)
  if(!is.null(list_extremes)){
    extrema = c(extrema, as.character(track_scales[list_extremes,ncol(track_scales)]))
  }
}
blacklist2 = unique(extrema)

all = all[!(all$trackID %in% blacklist2),]
alltest = all[,!(colnames(all)%in%remov.cols)]
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]
alltest[is.na(alltest)]=0
alltest = alltest[,- as.numeric(which(apply(alltest, 2, var) == 0))]

alltestscale = data.frame(scale(alltest))
#alltestscale = merge(cellfro,alltestscale,by=c("finalID","sample"))
alltestscale$type = all$type
alltestscale$frame = all$frame
alltestscale$last_frame = sapply(all$finalID,function(x) strsplit(x,"_")[[1]][1])


write.csv(x = alltestscale,file = "analysis/20210409_LGE_results.csv",quote = F, row.names = F, col.names = F)

#alltestscale = alltestscale[,- as.numeric(which(apply(alltestscale, 2, var) == 0))]
#corstuff = cor(alltestscale)
#heatmap(corstuff)


pca_alltest = prcomp((alltestscale))
plot((pca_alltest$sdev**2)/sum(pca_alltest$sdev**2),xlim = c(0,50))
plot(pca_alltest$x[,1],pca_alltest$x[,2])
pcaDF = data.frame(pca_alltest$x[,1:8])

#umap.defaults
umaptest = umap(alltestscale,metric='cosine')
um = data.frame(umaptest$layout); colnames(um) = c('U1','U2')
final = cbind(um,all,pcaDF)

ggplot(data=final[final$frame>0,],aes(U1, U2))+
  #geom_path(aes(fill=trackID))+
  #facet_wrap(~frame,ncol=4)+
  geom_point(aes(col=log2(CY5/YFP)),size=3)+
  scale_color_gradient2(low="blue",mid="grey",high="red",midpoint=0)+
# coord_cartesian(xlim=c(-20,20),ylim=c(-20,20))+
#  #geom_line(aes(fill=trackID,col=type))+
  theme_bw()

ggplot(final[final$frame>11,],aes(U1, U2))+
  facet_wrap(~sample,ncol=4)+
  geom_point(aes(col=log2(CY5norm/YFPnorm),size=magnitude))+
  geom_path(data= final[final$trackID %in% unique(final$trackID)[20],],aes(fill = trackID),size=2)+
  scale_color_gradient2(low="black",high="red")+
  #coord_cartesian(xlim=c(0,1050),ylim=c(0,1050))+
  #geom_line(aes(fill=trackID,col=type))+
  theme_bw()



ggplot(data=final,aes(frame,  area))+
  #facet_wrap(~sample,ncol=4)+
  #geom_point(aes(col=CY5,alpha=CY5))+
  geom_line(aes(fill=trackID,col=log2(CY5norm/YFPnorm)))+
  scale_color_gradient2(low="black",high="red")+
  theme_bw()

ggplot(data=final[final$frame>8,],aes(YFP, hara_7))+
  geom_point(aes(col=CY5,alpha=CY5))+
  #geom_line(aes(fill=trackID,col=CY5,alpha=CY5))+
  theme_bw()

late = final[final$frame>9,]

abcd = cor(data.matrix(all))

late = all[all$frame>9,]
early = all[all$frame<3,]


ggplot(data=all,aes(area))+
  facet_grid(type~.)+
  geom_histogram(aes(fill=type),col='darkgrey',binwidth=5)+
  theme_bw()

t.test(early[early$type == 'goblet',]$area,early[early$type == 'entero',]$area)

#####



featureIntensity = function(image.select){
  remove_str = gsub("\n","",image.select["intensity_image"])
  break_up = strsplit(remove_str, split = "]")
  tot.len = (length(break_up[[1]])-1)
  im.fin = c()
  for (i in 1:tot.len){
    remove_brk = gsub("\\[|\\]","",unlist(break_up)[[i]])
    num.vec = as.numeric(as.matrix(unlist(strsplit(gsub("\\s+", ",", gsub("^\\s+|\\s+$", "",remove_brk)), 
                                                   split = ","))))
    im.fin = rbind(im.fin,num.vec)
  }
  print(image(im.fin,col = gray.colors(12, rev = F)))
}

sampleGob = checking[107,]
CellID  = as.character(as.matrix(sampleGob$finalID))
experiment =  as.character(as.matrix(sampleGob$sample))

HCRDir = paste("results",experiment,"HCR","",sep = "/")
TL_HCR_connect = data.frame(read.csv(paste0(HCRDir,experiment,"_totalConnections.csv")))
TL_HCR_connect = TL_HCR_connect[TL_HCR_connect$Master_ID_13 == CellID,]
#featureIntensity(final[3,])
test= c()
if (nrow(TL_HCR_connect)>0){
  TL_characteristics  = as.matrix(TL_HCR_connect)[2:14]
  
  for(ii in 1:13){
    tmp = getStaticFeatures(TL_characteristics[ii],experiment = experiment, track =  as.character(TL_HCR_connect$TrackID), finalID=as.matrix(TL_HCR_connect)[ncol(TL_HCR_connect)])
    featureIntensity(tmp)
    test = rbind(test,tmp)
  }
}
