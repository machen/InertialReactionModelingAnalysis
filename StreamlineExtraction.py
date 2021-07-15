import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import re
import os
from mpl_toolkits.mplot3d import Axes3D  # Import to plot in 3D
import seaborn as sns  # Import to have nicer graphs

"""
extractStreamlines(data): Finds streamlines in data where the streamline
reverses in the y direction by establishing the initial direction
from the first point and taking the derivative.
UPDATE ME TO TAKE AN ARBITRARY GLOBAL DIRECTION.

Ideas:
-Develop reactive streamlines, determine fraction that is in recirculation pillar zone

To do: Read in and compare H2O2 streamlines with TCPO streamlines (perhaps we hash it or something to compare?)

1) Read in streamlines for the same file name as two separate data sets, determine using file extension
2) Streamline comparison-> does a given streamline ID correspond to the same set of coordinates to within some error?
3) Subselect streamline data to focus on pillar region
4) Select streamlines where Conc A AND Conc B exceeds 1% of initial
5) Discriminate selected streamlines between those experiencing any kind of recirculation, and those which do not
6) Calculate the fraction of streamlines which are in the vortex, plot vs the extracted Re number
7) Output: processing data is expensive, try to save yourself by outputting metadata file


"""


def subSelectData(data, xRange=None, yRange=None, zRange=None):
    """Assumes that each of the inputs to the function is a tuple containing
     max and min values of x, y, and z that we wish to include. Use for rough
     chopping of the data to exclude main channel flows"""
    if xRange:
        data = data.loc[(data.x > min(xRange)) & (data.x < max(xRange)), :]
    if yRange:
        data = data.loc[(data.y > min(yRange)) & (data.y < max(yRange)), :]
    if zRange:
        data = data.loc[(data.z > min(zRange)) & (data.z < max(zRange)), :]
    return data


def extractParams(fileName, nPil=2, caseExt='flowdata.txt'):
    # Produces a dictionary of experimental parameters
    rePat = re.compile('Re(\d+\.?\d*).')
    dPat = re.compile('d(\d+\.?\d*)_')
    cPat = re.compile('c(\d+\.?\d*)')
    kPat = re.compile('k(\d+\.?\d*)_')
    dVal = re.search(dPat, fileName).group(1)
    reVal = re.search(rePat, fileName).group(1)
    if caseExt == 'chemdata.txt':
        cVal = re.search(cPat, fileName).group(1)
        kVal = re.search(kPat, fileName).group(1)
    else:
        cVal = 0
        kVal = 0
    if nPil == 2:
        r1Pat = re.compile('r1_(\d+?)_')
        r2Pat = re.compile('r2_(\d+?)_')
        r1Val = re.search(r1Pat, fileName).group(1)
        r2Val = re.search(r2Pat, fileName).group(1)
        res = {'r1': float(r1Val), 'r2': float(r2Val), 'd': float(dVal),
               'Re': float(reVal), 'c': float(cVal), 'k': float(kVal)}
    if nPil == 1:
        rPat = re.compile('r(\d+?)_')
        rVal = re.search(rPat, fileName).group(1)
        res = {'r1': float(rVal), 'r2': float(rVal), 'd': float(dVal),
               'Re': float(reVal), 'c': float(cVal), 'k': float(kVal)}
    return res


def extractStreamlines(data, minLen=1):
    # data = data with x, y, z, sID, value
    # minLen = minimum length of streamline
    dData = data.diff()
    # Select for data where the y value reverses direction
    # Uses the first row to determine predominant flow direction
    if dData.loc[1, 'y'] > 0:
        subData = data.loc[dData.y < 0, :]
        # Filter out entries not on the same streamline
        subData = subData.loc[dData.sID == 0, :]
        # Pick entries that go beyond a certain Y value
    elif dData.loc[1, 'y'] < 0:
        subData = data.loc[dData.y > 0, :]
        subData = subData.loc[dData.sID == 0, :]
    # Create a slice of data containing only the streamlines targetted above
    # groups = subData.groupby(subData.sID)
    # tgt_sID = groups.sID.count()[groups.sID.count() > minLen].index.values
    tgt_sID = np.array(subData.sID.unique()) # Only selects unique values
    tgt_data = data.loc[data.sID.isin(tgt_sID), :]
    groups = tgt_data.groupby(tgt_data.sID)
    minsID = groups.sID.count()[groups.sID.count() > minLen].index.values
    print("No. Backward streamlines found: {}".format(len(minsID)))
    tgt_data = tgt_data.loc[tgt_data.sID.isin(minsID)]
    dtgt_data = dData.loc[data.sID.isin(minsID), :]
    # Captures the streamlines by using the transition
    startPoints = tgt_data.loc[dtgt_data.sID != 0, :]
    return tgt_sID, tgt_data, startPoints


def reactiveStreamlineExtraction(dataA, dataB, thresh,
                                 regionX=None, regionY=None, regionZ=None):
    """
    Outputs the streamline IDs which are reactive and the IDs which are within vortices
    """
    # Compare the streamlines in the two data sets
    # Subselect data to focus only on pillar region
    subDataA = subSelectData(dataA, regionX, regionY, regionZ)
    # Selects streamlines where both data sets meet the threshold value criteria
    subDataA['thresh'] = False
    subDataA.loc[subDataA.loc[:, 'val'] > thresh, 'thresh'] = True
    subDataB = subSelectData(dataB, regionX, regionY, regionZ)
    subDataB['thresh'] = False
    subDataB.loc[subDataA.loc[:, 'val'] > thresh, 'thresh'] = True
    reactiveID = subDataA.loc[(subDataA.val >= thresh) & (subDataB.val >= thresh), 'sID'].unique()
    dDataA = dataA.diff()
    dSubDataA = subDataA.diff()
    # Select for data where the y value reverses direction
    # Uses the first row to determine predominant flow direction
    if dDataA.loc[1, 'y'] > 0:
        subData = subDataA.loc[dSubDataA.y < 0, :]
        # Filter out entries not on the same streamline
        subData = subData.loc[dSubDataA.sID == 0, :]
        # Pick entries that go beyond a certain Y value
    elif dDataA.loc[1, 'y'] < 0:
        subData = subDataA.loc[dSubDataA.y > 0, :]
        subData = subData.loc[dSubDataA.sID == 0, :]
    vortIDA = np.array(subData.sID.unique())

    dDataB = dataB.diff()
    dSubDataB = subDataB.diff()
    if dDataB.loc[1, 'y'] > 0:
        subData = subDataB.loc[dSubDataB.y < 0, :]
        # Filter out entries not on the same streamline
        subData = subData.loc[dSubDataB.sID == 0, :]
        # Pick entries that go beyond a certain Y value
    elif dDataB.loc[1, 'y'] < 0:
        subData = subDataB.loc[dSubDataB.y > 0, :]
        subData = subData.loc[dSubDataB.sID == 0, :]
    vortIDB = np.array(subData.sID.unique())

    return reactiveID, vortIDA, vortIDB


def pillarGapCalculation(r1, r2, d):
    """ Use this to properly calculate the true "pillar gap" area depending on
    the available parameters of the model. This calculation will directly assume:
    1) That the edge of the most upstream pillar is at y = 0
    2) That there are only two pillars
    3) Positive y is upstream, negative y is downstream
    4) The simulation base unit is microns

    The gap is defined by the edges of the pillars along the channel length
    and the width of the largest pillar
    """
    if not r2:
        r2 = r1
    xPil = 250
    yPil = -(2*r1+d/2)
    x1 = xPil-max(r1, r2)
    x2 = xPil+max(r1, r2)
    y1 = yPil+d/2
    y2 = yPil-d/2

    return [x1, x2], [y1, y2]


# plt.ion()  # Keep interactive mode on
workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemStreamlines\\RawData\\"  # Directory you want to scan through
caseName = ""
extA = ".H2O2stream.txt"
extB = ".TCPOstream.txt"
plotData = True
saveData = False
exactGap = True

os.chdir(workingDir)
filePat = re.compile('(.*)('+extA+')')
fileList = os.listdir('.')
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2', 'd', 'Re', 'k', 'c'])
for f in fileList:
    if re.match(filePat, f):
        print(f)
        params = extractParams(f, nPil=1, caseExt='chemdata.txt')
        params['fileName'] = f
        # Need to think through loading in a file then finding its "partner"
        dataA = pd.read_table(f, sep='\s+', skiprows=8,
                              names=['x', 'y', 'z', 'sID', 'val'])
        fAlt = re.match(filePat, f).group(1)+extB
        dataB = pd.read_table(fAlt, sep='\s+', skiprows=8,
                              names=['x', 'y', 'z', 'sID', 'val'])
        # Need to ID params
        totalStreams = len(dataA.sID.unique())
        if exactGap:
            xRange, yRange = pillarGapCalculation(params['r1'], params['r2'], params['d'])
            outName = 'Exact Gap StreamlineMetaResults.csv'
        else:
            xRange = None
            yRange = None
            outName = 'Whole StreamlineMetaResults.csv'
        reactiveIDs, vortA, vortB = reactiveStreamlineExtraction(dataA, dataB, params['c']/100, regionX=xRange, regionY=yRange)
        """ # Okay, so we've gotten the IDs for reactive streamlines and the streamlines in vortices
        # We have to sub select the data that meets both.

        What am I outputting from this script?
        - Percentage of streamlines that are in vortex
        - Percentage of streamlines that are reactive
        - Percentage that are both?

        I also then want to plot them, perhaps in 2D?
        """
        vortData = dataA.loc[dataA.sID.isin(vortA), :]
        reacData = dataA.loc[dataA.sID.isin(reactiveIDs), :]
        reacVortData = reacData.loc[reacData.sID.isin(vortA), :]
        params['vortStream'] = len(vortData.sID.unique())
        params['reacStream'] = len(reacData.sID.unique())
        params['reacVortStream'] = len(reacVortData.sID.unique())
        params['totalStream'] = len(dataA.sID.unique())

        metaData = metaData.append(params, ignore_index=True)
        # if saveData:
        #     startPoints.to_csv(f[:-4]+"_startPoints.csv")
        if plotData:
            if params['Re'] != 100:
                continue
            if params['k'] != 2000:
                continue
        if plotData:
            fig = plt.figure()
            #ax = fig.add_subplot(111, projection='3d')
            ax = fig.add_subplot(111)
            plt.title(f)
            for ID in dataA.sID:
                # Plot the streamline in 3D space, color by streamline
                if ID not in reactiveIDs:
                    continue
                else:
                    if ID in vortA:
                        color = 'r'
                    else:
                        color = 'b'
                subData = dataA.loc[dataA.sID == ID, :]
                # ax.plot(dataA.loc[dataA.sID == ID, 'x'],
                #         dataA.loc[dataA.sID == ID, 'y'],
                #         dataA.loc[dataA.sID == ID, 'z'],
                #         color=color, ls='-')
                ax.plot(subData.x,
                        subData.y,
                        color=color, ls='-')
            ax.set_aspect('equal')
            ax.set_zlim([0, 100])
            ax.set_ylim([-550, -250])
            ax.set_xlim([150, 350])
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
if saveData:
    metaData.to_csv(outName)