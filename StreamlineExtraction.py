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

"""


def extractStreamlines(data):
    dData = data.diff()
    # Select for data where the y value reverses direction
    if dData.loc[1, 'y'] > 0:
        subData = data.loc[dData.y < 0, :]
    elif dData.loc[1, 'y'] < 0:
        subData = data.loc[dData.y > 0, :]
    # Filter out entries not on the same streamline
    subData = subData.loc[dData.sID == 0, :]

    # Create a slice of data containing only the streamlines targetted above
    tgt_sID = np.array(subData.sID.unique())
    print("No. Backward streamlines found: {}".format(len(tgt_sID)))
    tgt_data = data.loc[data.sID.isin(tgt_sID), :]
    dtgt_data = dData.loc[data.sID.isin(tgt_sID), :]
    startPoints = tgt_data.loc[dtgt_data.sID != 0, :]
    return tgt_sID, tgt_data, startPoints


plt.ion()  # Keep interactive mode on
workingDir = "..\\Two Pillar Studies\\CoarseResults\\"  # Directory you want to scan through
caseName = "TwoInletsTwoColumns_coarse"
tgtExt = ".txt"
plotData = True

# FILE READER VERSION, NO PLOTTING

os.chdir(workingDir)
filePat = re.compile('streamlines_'+caseName+'\w*'+tgtExt)
fileList = os.listdir('.')
for f in fileList:
    if re.match(filePat, f):
        print(f)
        data = pd.read_table(f, sep='\s+', skiprows=8,
                             names=['x', 'y', 'z', 'sID', 'val'])
        tgt_sID, tgt_data, startPoints = extractStreamlines(data)
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

# # plt.show()
# startPoints.to_csv(tgtFile+"_startPoints.csv")
