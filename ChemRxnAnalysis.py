import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns

""" Script is for processing chemical reaction data. I need to decide what I want
to do!

-Calculate total mass
-Cross sectional plotting or axial integration?
-Calculate scalar dissipation rate? -> wait actually do this

%TODO:

FORMAT IS:
x, y, z, MeshID, MeshVolumeScale, MeshElementVolume, u (m/s), v (m/s), w (m/s),
p, Velocity Magnitude (m/s), Mass flow (kg/s), c_H2O2 (mol/m^3),
c_TCPO (mol/m^3), c_product (mol/m^3), k (m^3/(s*mol))


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
        data = data.loc[(data.z > min(zRange)) & (data.z < min(zRange)), :]
    return data


def extractParams(fileName, nPil=2):
    # Produces a dictionary of experimental parameters
    rePat = re.compile('Re(.*?).txt')
    dPat = re.compile('d(\d+?)_')
    dVal = re.search(dPat, fileName).group(1)
    reVal = re.search(rePat, fileName).group(1)
    if nPil == 2:
        r1Pat = re.compile('r1_(\d+?)_')
        r2Pat = re.compile('r2_(\d+?)_')
        r1Val = re.search(r1Pat, fileName).group(1)
        r2Val = re.search(r2Pat, fileName).group(1)
        res = {'r1': float(r1Val), 'r2': float(r2Val), 'd': float(dVal),
               'Re': float(reVal)}
    if nPil == 1:
        rPat = re.compile('r(\d+?)_')
        rVal = re.search(rPat, fileName).group(1)
        res = {'r1': float(rVal), 'r2': float(rVal), 'd': float(dVal),
               'Re': float(reVal)}
    return res


# Read through files in a directory


workingDir = "..\\Comsol5.4\\Multipillar\\Normal\\FlowData_Normal\\"
# workingDir = "."
caseName = "Multipillar_v5.2_Normal"
caseExt = "\.txt$"

dataRegion = [-2500, 250]  # [-5000, 250]
nBins = 200
logBins = True  # True to use log spaced bins, False to use linear bins
nPil = 1  # Number of pillars in file specification

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')

# Check for existence of a metadata file which should be something like:
# caseName+"_meta.csv"
# Purpose of metadata file is to say what we've run already
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2',
                                     'd', 'Re', 'dP', 'q', 'l'])
for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        # Check for fileName already in metaData, skip if so
        data = pd.read_table(fileName, header=9, sep='\s+',
                             names=['x', 'y', 'z', 'meshID', 'EleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow'])
        data = subSelectData(data, yRange=dataRegion)

