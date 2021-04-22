import re
import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms.shortest_paths.weighted import single_source_dijkstra
from networkx.algorithms.simple_paths import all_simple_paths
import networkx.algorithms.shortest_paths.weighted as ssp


def TLGconvert(connectFilename):
    t0 = pd.read_csv(connectFilename) ;
    currentFrame = connectFilename.split("/")[-1].split(".")[0]
    nextFrame = str(int(currentFrame)+1)
    t0["frameNumber"] = int(nextFrame)+1
    t0["parentID"] = int(currentFrame)*cellBase + t0[currentFrame]
    t0["pointID"] = int(nextFrame)*cellBase + t0[nextFrame]
    return(t0[['pointID', 'parentID', 'cost', 'frameNumber']][:-1])

def generateSource(dataFrame,openingCost=25):
    subFrame = dataFrame[['frameNumber','parentID']].drop_duplicates()
    naughtFrame = pd.DataFrame()
    for i in range(0,len(subFrame)):
        tmpFrame = int(subFrame.iloc[i].frameNumber)
        tmpPoint = int(subFrame.iloc[i].parentID)
        tmpCost = int(tmpFrame-2)*openingCost+1
        tmpFrame = pd.DataFrame(np.array([[tmpPoint,0,tmpCost,tmpFrame-1]]),
                                columns=['pointID', 'parentID', 'cost', 'frameNumber'])
        naughtFrame = naughtFrame.append(tmpFrame)
    return(naughtFrame)

def generateSink(dataFrame,closingCost=25):
    subFrame = dataFrame[['frameNumber','pointID']].drop_duplicates()
    finalFrame = subFrame['frameNumber'].max()+1
    finalPoint = int(10**np.ceil(np.log10(subFrame['pointID'].max())))
    omegaFrame = pd.DataFrame()

    for i in range(0,len(subFrame)):
        tmpFrame = int(subFrame.iloc[i].frameNumber)
        tmpPoint = subFrame.iloc[i].pointID
        tmpCost = (finalFrame-tmpFrame-1)*closingCost+1
        tmpFrame = pd.DataFrame(np.array([[finalPoint,tmpPoint,tmpCost,finalFrame]]),
                                    columns=['pointID', 'parentID', 'cost', 'frameNumber'])
        omegaFrame = omegaFrame.append(tmpFrame)
    return(omegaFrame)

# generate node dataFrame
def pointCapacity(dataFrame, divTime):
    frameA = dataFrame[['pointID','frameNumber']].drop_duplicates()
    frameA.loc[:,'frameNumber'] = frameA.loc[:,'frameNumber']+1
    frameV = dataFrame[['parentID','frameNumber']].drop_duplicates()
    frameV = frameV.rename(columns={"parentID": "pointID"})
    pointCapacities = frameV.append(frameA).drop_duplicates()
    pointCapacities['capacity'] = (np.floor(((pointCapacities['frameNumber'].max()-pointCapacities['frameNumber'])+(divTime-1))/divTime)+1).astype(int)
    pointCapacities['usage'] = 0
    return(pointCapacities)

def shortestSuccessivePath(originalFrame, openingCost = 100, closingCost = 100, maxCost=105):

    ptCapFrame = pointCapacity(originalFrame, divTime=10)

    openingFrame = generateSource(originalFrame,openingCost=openingCost)
    closingFrame = generateSink(originalFrame,closingCost=closingCost)
    newFrame = openingFrame.append(originalFrame.append(closingFrame))
    newFrame = newFrame.reset_index(drop=True)
    newFrame['weight'] = newFrame.cost

    newFrame['capacity'] = (newFrame['frameNumber'].max()-newFrame.frameNumber)
    newFrame.loc[newFrame['capacity']>1,'capacity']=2
    newFrame.loc[newFrame['capacity']<=1,'capacity']=1


    newFrame['usage'] = 0
    G=nx.from_pandas_edgelist(newFrame, 'parentID', 'pointID', ['weight','capacity','usage'],create_using=nx.DiGraph)
    maxPt = newFrame['pointID'].max()
    totalCost = 0
    ticker=0
    while totalCost < maxCost:
        tmp_shortest_cost = ssp.bellman_ford_path_length(G, 0, maxPt)
        totalCost = tmp_shortest_cost
        tmp_shortest_path = ssp.bellman_ford_path(G, 0, maxPt)
        for i in range(0,len(tmp_shortest_path)-2):
            #print( tmp_shortest_path[i+1])
            #ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i],'usage'] =( ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i],'usage']+ 1)
            ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i+1],'usage'] =( ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i+1],'usage']+ 1)

            #truthState0 = np.asarray(ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i],'usage'] >= ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i],'capacity'])[0]
            truthState1 = np.asarray(ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i+1],'usage'] >= ptCapFrame.loc[ptCapFrame.pointID == tmp_shortest_path[i+1],'capacity'])[0]
            #if truthState0:
            #    G.remove_node(tmp_shortest_path[i])
            if truthState1:
                G.remove_node(tmp_shortest_path[i+1])

        ticker=ticker+1
        print(ticker, ":",round(tmp_shortest_cost),":", tmp_shortest_path)

    return(ptCapFrame)
