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


def extractParams(fileName, nPil=2, caseExt='flowdata.txt'):
    # Produces a dictionary of experimental parameters
    rePat = re.compile('Re(\d+\.?\d*).'+caseExt)
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
    # data = data with x, y, z, sID
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


plt.ion()  # Keep interactive mode on
workingDir = "..\\Comsol5.3\\Two Pillar Studies\\CoarseResults\\v1Results\\"  # Directory you want to scan through
caseName = "TwoInletsTwoColumns_coarse"
tgtExt = ".txt"
extA = "H2O2stream.txt"
extB = "TCPOstream.txt"
plotData = True
saveData = False

# FILE READER VERSION, NO PLOTTING

os.chdir(workingDir)
filePat = re.compile('streamlines_'+caseName+'\w*'+tgtExt)
fileList = os.listdir('.')
for f in fileList:
    if re.match(filePat, f):
        print(f)
        data = pd.read_table(f, sep='\s+', skiprows=8,
                             names=['x', 'y', 'z', 'sID', 'val'])
        tgt_sID, tgt_data, startPoints = extractStreamlines(data, 100)
        if saveData:
            startPoints.to_csv(f[:-4]+"_startPoints.csv")
        if plotData:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            plt.title(f)
            for ID in tgt_sID:
                # Plot the streamline in 3D space, color by streamline
                subData = tgt_data.loc[tgt_data.sID == ID, :]
                ax.plot(subData.x, subData.y, subData.z, label=ID)
                ax.set_zlim([0, 100])
                ax.set_ylim([-5000, 2000])
                ax.set_xlim([-2000, 2000])
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                ax.set_zlabel('z')
