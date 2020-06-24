import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os

""" Purpose of script is to import data for velocity analysis,
but will also serve as a means of doing general data import from comsol.

%TODO:

-Can I natively capture the column names? -> (\w+? \(.+?\)|\w+?) @ \d+?: will
capture the header name, but it will not capture the resultant parameter info.
Likely that info is not worth it, and should be kept in the file name.
BUT SHOULD I? No. I should format my data consistently

FORMAT IS:
x, y, z, MeshID, MeshVolumeScale, MeshElementVolume, u (m/s), v (m/s), w (m/s),
p, Velocity Magnitude (m/s), Mass flow (kg/s)

-Analyze for recirculation zone (how do we tag points as in or not in a zone?)
-Do calculations for residence time?
    How?

-Metadata file handling:

-Export
"""


def produceVelPDF(data, nBins):
    """OUTPUT"
    normFreq: Normalized frequencies for the defined bins
    velVal: The mean velocity in each bin
    groups: The binned data group
    velBin: Velocity bin limits used
    Produce a velocity histogram from the volumes and velocities
    1) Bin the velocities by magnitude
    2) Find the total grid size in each bin
    3) Weight the frequency by multiplying the number of bin entries
    by the total vol in the bin
    4) Normalize the frequency to make the area under the curve 1.
    (I think it's sum (freq*area))
    # Bin velocities -> is there a better way to do that?
    We might need it since the bins need to be better at low velocities
    # Option, bin velocities such that the mean and median are within
    some tolerance or the RSD meets a certain criteria
    # velBin = np.linspace(data.velMag.min(), data.velMag.max(), num=nBins)"""
    velBin = np.logspace(np.log10(data.velMag.min()),
                         np.log10(data.velMag.max()), num=nBins)
    groups = data.groupby(np.digitize(data.velMag, velBin))
    velVal = groups.velMag.mean()
    # Weight frequencies
    weightedFreq = groups.EleVol.sum()*groups.size()
    # Calculate area under current curve
    totalArea = np.trapz(weightedFreq)
    normFreq = weightedFreq/totalArea
    return normFreq, velVal, groups, velBin


def extractParams(fileName):
    # Produces a dictionary of experimental parameters
    r1Pat = re.compile('r1_(\d+?)_')
    r2Pat = re.compile('r2_(\d+?)_')
    rePat = re.compile('Re(.*?).txt')
    dPat = re.compile('d(\d+?)_')
    r1Val = re.search(r1Pat, fileName).group(1)
    r2Val = re.search(r2Pat, fileName).group(1)
    dVal = re.search(dPat, fileName).group(1)
    reVal = re.search(rePat, fileName).group(1)
    res = {'r1': float(r1Val), 'r2': float(r2Val), 'd': float(dVal),
           'Re': float(reVal)}
    return res


def calcFlowPress(data, params, nu=1.6E-6, c=500E-6, cRatio=0.5,
                  depth=100E-6):
    # vel = Re*nu/(c*cRatio)
    # q0 = vel*c*depth
    # nu is default kinematic viscosity for acetonitrile
    # c is the channel width
    # cRatio is the ratio of the inlet channel width to the outlet channel
    # depth is the channel depth
    velInlet = params['Re']*nu/(c*cRatio)
    q0 = velInlet*c*depth  # m^3/s
    deltaP = data.p.max()-data.p.min()  # Pa
    return deltaP, q0


# Read through files in a directory


workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\Normal\\FlowData\\"
caseName = "TwoInletsTwoColumns_v5.1_Normal_FlowData"
caseExt = ".txt"
createFresh = True  # Create a new metadata file

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')
# Check for existence of a metadata file which should be something like:
# caseName+"_meta.csv"
# Purpose of metadata file is to say what we've run already
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2',
                                     'd', 'Re', 'dP', 'q'])
for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        # Check for fileName already in metaData, skip if so
        data = pd.read_table(fileName, header=9, sep='\s+',
                             names=['x', 'y', 'z', 'meshID', 'EleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow'])
        # There's no need to do any of this because produceVelPDF does it
        # avgVol = data['EleVol'].mean()
        # data['NormScale'] = data['EleVol']/avgVol
        # # The velMag doesn't need to be calculated, since it's output
        # # data['velMagCalc'] = np.sqrt(data.u**2 + data.v**2 + data.w**2)
        # # data['velMagErr'] = (data.velMagCalc-data.velMag)/data.velMag
        # data['velMagScaled'] = data.velMag/data.NormScale
        # data['uScaled'] = data.u/data.NormScale
        # data['vScaled'] = data.v/data.NormScale
        # data['wScaled'] = data.w.values/data.NormScale.values
        normFreq, velVals, velGroups, velBin = produceVelPDF(data, 1000)
        params = extractParams(fileName)
        params['dP'], params['q'] = calcFlowPress(data, params)
        params['fileName'] = fileName
        metaData = metaData.append(params, ignore_index=True)
        velData = {'NormFreq': normFreq, 'velVal': velVals}
        velPDF = pd.DataFrame(velData)
        velPDF.to_csv(fileName[:-4]+"_histogram.csv")
        plt.figure()
        plt.plot(velVals, normFreq)
        plt.xlabel('Average velocity of bin (m/s)')
        plt.ylabel('Normalized Frequency (.)')
        plt.savefig(fileName[:-4]+".png")
        plt.close()

metaData.to_csv(caseName+"_meta.csv")

# Plot deltaP vs Q, also fit a line
plt.plot(metaData.dP, metaData.q, ls='None', marker='*')
linFit = np.polyfit(metaData.dP, metaData.q, 1)
interpQ = np.polyval(linFit, metaData.dP)
plt.plot(metaData.dP, interpQ, ls='-', color='k')
plt.xlabel('Pressure Difference (Pa)')
plt.ylabel('Flow rate (m^3/s)')
plt.title('Pressure vs Flow Rate')
plt.savefig(caseName+'_dPvsQ.png')
