from genericpath import isdir
from readline import set_completion_display_matches_hook
import DataHelper as dh
import os, re
import pandas as pd
"""Script is meant to do calculate values of things using integrated quantities
(vol. weighted conc., pressure, etc.) in slices, and outputting those results
to a readable spreadsheet"""


def genSliceFromPillar(data, xPil, yPil, rPil, split=True):
    # Should generate coordinates that can be used to get slices for an XZ plane
    # Slice can be around the pillar, or inlcude the
    zMax = data.z.max()
    zMin = data.z.min()
    xMax = data.x.max()
    xMin = data.x.min()
    if not split:
        # Assumes you want the whole sim just sliced at one line
        slice = [xMin, xMax, zMin, zMax, yPil]
        return slice
    else:
        # Slices above and below the pillar, only works for a single line of pillars
        topSlice = [xPil+rPil, xMax, zMin, zMax, yPil]
        botSlice = [xMin, xPil-rPil, zMin, zMax, yPil]
    return topSlice, botSlice


dataDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\RawData\\"
workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\SliceData\\"
caseName = 'TwoPillar_v6'
caseExt = r'chemdata\.txt'
nPil = 2
sliceList = True


filePat = re.compile(caseName+'.*?'+caseExt+'$')
if not os.isdir(workingDir):
    raise FileNotFoundError('Please create working directory with slice specification')
if sliceList:
    # load slice list from directory. Useful for cases with the same geometry
    try:
        slices = pd.read_csv(workingDir+'Slices.csv', header=0)
    except FileNotFoundError:
        print('No slice specification, stopping script.')
        quit()


dataList = os.listdir(dataDir)
workingList = os.listdir(workingDir)

# Iterate through files, picking the ones matching the caseName and extension
for fileName in dataList:
    if re.match(filePat, fileName):
        print(fileName)
        data = dh.dataLoader(fileName, caseExt)
        params = dh.extractParams(fileName, nPil, caseExt)
        if not sliceList:
            slices = genSlices()
            # Create slice list from geometry of pillars. Need helper function to calculate pillar locations from parameters.