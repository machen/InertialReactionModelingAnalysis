import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns

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


def produceVelPDF(data, nBins=1000, logBin=True, prop="velMag"):
    """
    INPUT:
    data: data output from comsol, must include column specified by prop & "eleVol" cols
    nBins: Number of bins you want to use, default is 1000
    logBin: True-> bins are evenly spaced in log space,
    otherwise, linear spacing
    prop: Column name of data you wish to use, default, velMag (velocity magnitude)

    OUTPUT
    normFreq: Normalized frequencies for the defined bins by bin size
    velVal: The mean velocity in each bin
    groups: The binned data group
    velBin: Velocity bin limits used
    Produce a velocity PDF from the input data

    Option, bin velocities such that the mean and median are within
    some tolerance or the RSD meets a certain criteria

    Max velocity is taken as 1.01*max due to the way np.digitize works.
    For a linear bin, without doing this, the last value will not belong to a
    bin that has a proper size definition.
    """
    if logBin:
        velBin = np.logspace(np.log10(data.loc[:, prop].min()),
                             np.log10(data.loc[:, prop].max()*1.01), num=nBins)
    else:
        velBin = np.linspace(data.loc[:, prop].min(),
                             data.loc[:, prop].max()*1.01, num=nBins)
    velBinSize = velBin[1:]-velBin[:-1]
    data.loc[:, 'binID'] = np.digitize(data.loc[:, prop], velBin)
    groups = data.groupby(data.binID)
    velVal = groups[prop].mean()
    # Weight frequencies by included volume and normalize to bin size
    weightedFreq = groups.EleVol.sum()*groups.size() \
        / groups.binID.median().apply(lambda x: velBinSize[x-1])
    # Calculate area under freq-vel curve to obtain a PDF
    totalArea = np.trapz(weightedFreq, x=velVal)
    normFreq = weightedFreq/totalArea
    return normFreq, velVal, groups, velBin


def extractParams(fileName, nPil = 2):
    # Produces a dictionary of experimental parameters
    rePat = re.compile('Re(.*?).flowdata.txt')
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


def calcFlowPress(data, params, nu=1.6E-6, c=500E-6, cRatio=0.5,
                  depth=100E-6, regionSize=10):
    """# vel = Re*nu/(c*cRatio)
    # q0 = vel*c*depth
    # nu is default kinematic viscosity for acetonitrile
    # c is the channel width
    # cRatio is the ratio of the inlet channel width to the outlet channel
    # depth is the channel depth
    # Region size is how wide a slice you want to use for pressure averaging
    You should more explicitly calculate pressure of the inflow and outflow
    1) Group data by y value (could be 1 micron intervals) -> data.groupby(data.y)
    Or you could do subDataMin = data.loc[(data.y <= TARGETMAX) & (data.y>= TARGETMIN), :]
    2) find average pressure in each bin which represents a y slice-> WORTH OUTPUTTING POTENTIALLY
    3) Find absolute value of pressure difference of max y bin and min y bin
    4) Output some statistics about that pressure/info about where it's from
    maybe also the distance over which the pressure drop is realized"""
    velInlet = params['Re']*nu/(c*cRatio)
    q0 = velInlet*c*depth  # m^3/s
    yInlet = data.y.min()
    inletData = data.loc[(data.y <= yInlet+regionSize), :]
    yOutlet = data.y.max()
    outletData = data.loc[(data.y >= yOutlet-regionSize)]
    inletP = inletData.p.mean()
    outletP = outletData.p.mean()
    deltaP = abs(inletP - outletP)  # Pa
    dx = abs(inletData.y.mean() - outletData.y.mean())
    return deltaP, q0, dx


# Read through files in a directory


workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData_FlowOnly\\Raw Data\\"
# workingDir = "."
caseName = "TwoInletsTwoColumns_v5."
caseExt = "\.flowdata.txt$"
writeMeta = True  # Create new metadata files
binVel = True  # True to bin velocties, false to skip
binProp = 'v'  # Name of column to run PDF on

dataRegionX = [100, 400]
dataRegionY = [-700, -100]  # [-5000, 250]
nBins = 100
logBins = True  # True to use log spaced bins, False to use linear bins
nPil = 2  # Number of pillars in file specification

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
        data = subSelectData(data, xRange=dataRegionX, yRange=dataRegionY)
        params = extractParams(fileName, nPil)
        params['dP'], params['q'], params['l'] = calcFlowPress(data, params)
        params['fileName'] = fileName
        metaData = metaData.append(params, ignore_index=True)
        if binVel:
            normFreq, velVals, velGroups, velBin = \
                produceVelPDF(data, nBins=nBins, logBin=logBins, prop=binProp)
            velData = {'normFreq': normFreq, 'velVal': velVals}
            velPDF = pd.DataFrame(velData)
            velPDF.to_csv(fileName[:-4]+"_histogram.csv")
            plt.figure()
            plt.plot(velVals, normFreq)
            plt.xlabel('Average velocity of bin (m/s)')
            plt.ylabel('Normalized Frequency (.)')
            plt.savefig(fileName[:-4]+"_linear.png")
            plt.yscale('log')
            plt.xscale('log')
            plt.savefig(fileName[:-4]+"_log.png")
            plt.close()

flowFitData = pd.DataFrame([], columns=['r1', 'r2', 'd', 'linA', 'linB',
                                        'quadA', 'quadB', 'quadC', 'expA',
                                        'expB'])
# Obtain the unique geometries used (i.e. by d, r1, and r2)
uniqueParam = metaData.loc[~metaData.duplicated(['d', 'r1', 'r2']), :]
colorPalette = sns.color_palette('deep', n_colors=len(uniqueParam.index))

f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
ci = 0
for i in uniqueParam.index:
    r1 = uniqueParam.loc[i, 'r1']
    r2 = uniqueParam.loc[i, 'r2']
    d = uniqueParam.loc[i, 'd']
    subData = metaData.loc[(metaData.d == d) &
                           (metaData.r1 == r1) & (metaData.r2 == r2), :]
    ax1.plot(subData.q, subData.dP, ls='None', color=colorPalette[ci],
             marker='o', label='r1 = {}, r2 = {}, d = {}'.format(r1, r2, d))
    linFit = np.polyfit(subData.q, subData.dP, 1)
    quadFit = np.polyfit(subData.q, subData.dP, 2)
    logFit = np.polyfit(np.log(subData.q), np.log(subData.dP), 1)
    interpQ = np.linspace(min(subData.q), max(subData.q), num=100)
    interp_dP = np.polyval(linFit, interpQ)
    interpQuad_dP = np.polyval(quadFit, interpQ)
    interpLog_dP = np.exp(np.polyval(logFit, np.log(interpQ)))
    caseParam = {'r1': r1, 'r2': r2, 'd': d, 'linA': linFit[0],
                 'linB': linFit[1], 'quadA': quadFit[0], 'quadB': quadFit[1],
                 'quadC': quadFit[2], 'expA': logFit[0], 'expB': logFit[1]}
    flowFitData = flowFitData.append(caseParam, ignore_index=True)
    ax1.plot(interpQ, interp_dP, ls='-', color=colorPalette[ci], label="{:.2e}*q+{:2e}".format(*linFit))
    ax1.plot(interpQ, interpQuad_dP, ls='--', color=colorPalette[ci], label="{:.2e}*q^2+{:.2e}*q+{:.2e}".format(*quadFit))
    ax1.plot(interpQ, interpLog_dP, ls='dotted', label="log(dP) = {:.2e}*log(q)+{:.2e}".format(*logFit))
    ci += 1

if writeMeta:
    metaData.to_csv(caseName+"_meta.csv")
    flowFitData.to_csv(caseName+"_flowFits.csv")

"""
What do I want to do?

For each unique combination of d, r1, and r2, I want to apply a coloration
and marker to differentiate it, fit Q vs dP, and plot both.

Option 1: Create table of the unique entries
subselect data for each known entry

Option 2: Iterate through d, r1, r2 unique values
Bad because I may not have all combinations.
"""

# Plot deltaP vs Q, also fit a line
# plt.plot(metaData.q, metaData.dP, ls='None', marker='*', color='k',
#          label='Overall')
# linFit = np.polyfit(metaData.q, metaData.dP, 1)
# interp_dP = np.polyval(linFit, metaData.q)
# plt.plot(metaData.q, interp_dP, ls='-', color='k', label='Overall fit')
plt.legend(loc=0)
plt.ylabel('Pressure Difference (Pa)')
plt.xlabel('Flow rate (m^3/s)')
plt.title('Flow rate vs Pressure')
plt.savefig(caseName+'_dPvsQ.png', dpi=600)
