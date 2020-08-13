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


def dataReduction(data, xBin, yBin, zBin):
    """Reduce the data by averaging into a coraser grid.
    Should maybe also weight the entries by volume.
    How do we do this? I bet we can use the group by function to group things into
    a defined grid (like I want 50 entries in x y and z)"""

    return data


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


workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\ChemData\\"
# workingDir = "."
caseName = "TwoInletsTwoColumns_v5.2_ExF"
caseExt = "\.chemdata\.txt$"
writeMeta = True
dataRegionY = [-250, 1500]

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')

# Check for existence of a metadata file which should be something like:
# caseName+"_meta.csv"
# Purpose of metadata file is to say what we've run already
metaData = pd.DataFrame([], columns=['fileName', 'totalProd', 'totalH2O2',
                                     'totalTCPO', 'dCdtAvg', 'totalVol',
                                     'M'])
for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        # Check for fileName already in metaData, skip if so
        data = pd.read_table(fileName, header=10, sep='\s+',
                             names=['x', 'y', 'z', 'meshID', 'eleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow',
                                    'h2o2', 'tcpo', 'cProduct', 'k'])
        data = subSelectData(data, yRange=dataRegionY)
        data['dCdt'] = data.h2o2*data.tcpo*data.k
        data['prodMass'] = data.cProduct*data.eleVol
        data['tcpoMass'] = data.tcpo*data.eleVol
        data['h2o2Mass'] = data.h2o2*data.eleVol
        data['M'] = (data.tcpo+data.cProduct)*(data.tcpo+data.cProduct)*data.eleVol
        res = {}
        res['fileName'] = fileName
        res['totalProd'] = data.prodMass.sum()
        res['totalH2O2'] = data.h2o2Mass.sum()
        res['totalTCPO'] = data.tcpoMass.sum()
        res['totalVol'] = data.eleVol.sum()
        res['dCdtAvg'] = data.dCdt.mean()
        res['M'] = data.M.sum()

        # data = subSelectData(data, yRange=dataRegionY)
        metaData = metaData.append(res, ignore_index=True)

if writeMeta:
    metaData.to_csv(caseName+"_meta.csv")
