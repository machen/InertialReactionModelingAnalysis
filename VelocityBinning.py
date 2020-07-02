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


def produceVelPDF(data, nBins=1000, logBin=True):
    """
    INPUT:
    data: data output from comsol, must include "velMag" & "eleVol" cols
    nBins: Number of bins you want to use, default is 1000
    logBin: True-> bins are evenly spaced in log space,
    otherwise, linear spacing

    OUTPUT
    normFreq: Normalized frequencies for the defined bins by bin size
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

    Option, bin velocities such that the mean and median are within
    some tolerance or the RSD meets a certain criteria

    Bin size is taken as 1.01*max due to the way np.digitize works.
    For a linear bin, without doing this, the last value will not belong to a
    bin that has a proper size definition.
    """
    if logBin:
        velBin = np.logspace(np.log10(data.velMag.min()),
                             np.log10(data.velMag.max()*1.01), num=nBins)
    else:
        velBin = np.linspace(data.velMag.min(), data.velMag.max()*1.01,
                             num=nBins)
    velBinSize = velBin[1:]-velBin[:-1]
    data.loc[:, 'binID'] = np.digitize(data.velMag, velBin)
    groups = data.groupby(data.binID)
    velVal = groups.velMag.mean()
    # Weight frequencies by included volume and normalize to bin size
    weightedFreq = groups.EleVol.sum()*groups.size() \
        / groups.binID.median().apply(lambda x: velBinSize[x-1])
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


# workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData_FlowOnly\\"
workingDir = "."
caseName = "TwoInletsTwoColumns_v5.1_Normal_FlowData"
caseExt = "\.txt$"
writeMeta = True  # Create a new metadata file
binVel = True  # True to bin velocties, false to skip

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
        params = extractParams(fileName)
        params['dP'], params['q'] = calcFlowPress(data, params)
        params['fileName'] = fileName
        metaData = metaData.append(params, ignore_index=True)
        if binVel:
            normFreq, velVals, velGroups, velBin = \
                produceVelPDF(data, nBins=1000, logBin=True)
            velData = {'NormFreq': normFreq, 'velVal': velVals}
            velPDF = pd.DataFrame(velData)
            velPDF.to_csv(fileName[:-4]+"_histogram.csv")
            plt.figure()
            plt.plot(velVals, normFreq)
            plt.xlabel('Average velocity of bin (m/s)')
            plt.ylabel('Normalized Frequency (.)')
            plt.savefig(fileName[:-4]+".png")
            plt.close()

if writeMeta:
    metaData.to_csv(caseName+"_meta.csv")

# Obtain the unique geometries used (i.e. by d, r1, and r2)
# metaData.loc[~metaData.duplicated(['d', 'r1', 'r2']), :]

# Plot deltaP vs Q, also fit a line
plt.plot(metaData.q, metaData.dP, ls='None', marker='*')
linFit = np.polyfit(metaData.q, metaData.dP, 1)
interp_dP = np.polyval(linFit, metaData.q)
plt.plot(metaData.q, interp_dP, ls='-', color='k')
plt.xlabel('Pressure Difference (Pa)')
plt.ylabel('Flow rate (m^3/s)')
plt.title('Pressure vs Flow Rate')
plt.savefig(caseName+'_dPvsQ.png')
