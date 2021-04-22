import re
import pandas as pd
import numpy as np
import argparse


def tlgHistoryRecall(table, pointID):
    init = table[table.pointID == pointID]
    track = ":"+str(init.frameNumber.values[0])+"-"+str(init.pointID.values[0])
    parental = init.parentID.values[0]
    while (parental>0) :
        init = table[table.pointID == parental]
        track = ","+str(init.frameNumber.values[0])+"-"+str(init.pointID.values[0]) + track
        parental = init.parentID.values[0]

    return(track)

def getAncestory(table,pointID):
    init = table[table.pointID == pointID]
    parental = init.parentID.values[0]
    div=init.Div.values[0]
    ttLast = 0
    if (div==1):
        div=0;
    while (parental>0)&(div ==0) :
        init = table[table.pointID == parental]
        parental = init.parentID.values[0]
        div=init.Div.values[0]
        ttLast = ttLast+1
    lcAncestor = init.pointID.values[0]
    outFrame = pd.DataFrame({"pointID": [pointID], "timefromDiv": [ttLast], "lcAncest": [lcAncestor]})
    return(outFrame)

def detailTLGoutput(table):
    parents =  table.parentID.unique()[table.parentID.unique()!=0]
    addon1 = pd.DataFrame()
    for i in parents:
        checklen = (len(table[table.parentID == i]))
        if(checklen==1):
            div = 0
        else:
            div = 1
        addon1 = addon1.append(pd.DataFrame({"pointID": [i],"Div": [div]}))

    points = table.pointID.unique()
    res = [i for i in points if i not in parents]
    addon2 = pd.DataFrame({"pointID": res, "end": 1})

    tst2 = pd.merge(pd.merge(table,addon1,how="outer"),addon2,how="outer")
    tst2 = tst2.fillna(0)
    tst2 = tst2.astype({'Div': 'int', 'end': 'int'})

    nde = pd.DataFrame()
    for j in tst2.pointID.values:
        tmp_ancestor = getAncestory(tst2,j)
        nde = nde.append(tmp_ancestor)

    final = pd.merge(tst2,nde,how="outer")
    return(final)

def ancestorStr(refTable,ancestor):
    tmp=refTable[refTable.lcAncest==ancestor];
    if len(tmp.pointID.values)>1:
        tmp = tmp[tmp.lcAncest != tmp.pointID]
        if len(tmp.pointID.values)>1:
            str0 = str(tmp.pointID.values[0])+":"+str(tmp.timefromDiv.values[0])
            str1 = str(tmp.pointID.values[1])+":"+str(tmp.timefromDiv.values[1])
            return("("+str0+","+str1+")"+str(ancestor))
        else:
            str0 = str(tmp.pointID.values[0])+":"+str(tmp.timefromDiv.values[0])
            return("("+str0+")"+str(ancestor))
    else:
        str0 = str(tmp.pointID.values[0])+":"+str(tmp.timefromDiv.values[0])
        return("("+str0+")"+str(ancestor))

def generateNewick(TLGfile):

    # Load output from TimelapseGUI
    tstng = pd.read_csv(TLGfile) ;

    # Use function to generate extra data
    tlg = detailTLGoutput(tstng)

    # Determine the ancestor and end nodes
    ancestorNodes = tlg.lcAncest.unique();
    endNodes = tlg[tlg.end==1].pointID.unique()

    # Refine the search to only include ancestors and ends
    refined = tlg[tlg['pointID'].isin(np.concatenate([ancestorNodes,endNodes]))]

    refined = refined.fillna(0)
    refined = refined.astype({'lcAncest': 'int', 'pointID': 'int','parentID': 'int','frameNumber': 'int'})
    refined.frameNumber = refined.frameNumber+1
    #refined.loc[refined.parentID == 0,'Div'] = 1
    #refined.loc[refined.parentID == 0,'lcAncest'] = 0
    #refined = refined[:-1]

    refined['coupling'] = refined['lcAncest'].apply(lambda x: ancestorStr(refined,x))

    toplvlIDs=refined[(refined.lcAncest == refined.pointID)&(refined.parentID == 0)][:-1].pointID.values
    replaceForTop = "("+":1,".join(str(x) for x in toplvlIDs)+":1)0;"
    refined.loc[refined.lcAncest==refined.pointID,'coupling'] = replaceForTop

    refined = refined.loc[~refined['pointID'].isin(toplvlIDs)]

    LCAdf = pd.DataFrame() ;

    for i in refined.lcAncest.unique():
        tmpLCA = refined[refined.lcAncest == i]
        tmpLCA = tmpLCA.loc[:,['lcAncest','frameNumber','coupling']].drop_duplicates()
        LCAdf = LCAdf.append(tmpLCA[tmpLCA.frameNumber == tmpLCA.frameNumber.min()])

    LCAsort = LCAdf.sort_values(by=['frameNumber'])
    initial_line = LCAsort.iloc[0].coupling

    for i in range(1,len(LCAsort)):
        pattern = "\("+str(LCAsort.iloc[i].lcAncest)+":"
        replace = "("+LCAsort.iloc[i].coupling+":"
        newline = re.sub(pattern, replace,initial_line)
        if(newline == initial_line):
            pattern = "\,"+str(LCAsort.iloc[i].lcAncest)+":"
            replace = ","+LCAsort.iloc[i].coupling+":"
            newline = re.sub(pattern, replace,initial_line)
        initial_line = newline

    return(initial_line)



def main():
    parser = argparse.ArgumentParser('Generate Newick-formatter tree.')

    parser.add_argument('--in', type = str, dest='file', help = 'Path and file in timelapse GUI output form (CSV).')

    args = parser.parse_args()

    newick = generateNewick(args.file)

    print(newick)

    return

if __name__ == '__main__':
    main()
