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
-This should report the left and right bin boundaries,
rather than doing interpolation of the values for the velocity of a given bin
Or one bin side and the corresponding bin size.


FORMAT IS:
x, y, z, MeshID, MeshVolumeScale, MeshElementVolume, u (m/s), v (m/s), w (m/s),
p, Velocity Magnitude (m/s), Mass flow (kg/s)

-Analyze for recirculation zone (how do we tag points as in or not in a zone?)
-Do calculations for residence time?
    How?

-Metadata file handling:

-Export
"""


def dataLoader(filename, type='flowdata.txt'):
    if type == 'flowdata.txt':
        data = pd.read_table(fileName, header=9, sep='\s+',
                             names=['x', 'y', 'z', 'meshID', 'EleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow',
                                    'vortX', 'vortY', 'vortZ', 'vortMag'])
    if type == 'chemdata.txt':
        data = pd.read_table(fileName, header=10, sep='\s+',
                             names=['x', 'y', 'z', 'meshID', 'eleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow',
                                    'h2o2', 'tcpo', 'cProduct', 'k'])
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


def crossings_nonzero_all(data):
    # Will calculate a zero crossing and return the indices where the crossing occurs
    # Credit to https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
    pos = data > 0
    npos = ~pos
    return ((pos[:-1] & npos[1:]) | (npos[:-1] & pos[1:])).nonzero()[0]


def producePDF(data, nBins=1000, logBin=True, prop="velMag"):
    """
    TO DO: Try to boost computational efficiency by subselecting the data to work on
    (i.e., don't use the columns you don't need)
    subData = data.loc[:,('EleVol', prop)]

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

    Alternate bin scheme to enable log binning of negative numbers and ranges
    that cross 0:
    -ID minimum absolute value and maximum absolute value
    -logbin that min to max range
    -duplicate range, reverse it, and then multiply by -1
    -concatenate the two ranges, including 0 between the two
    -pass the ranges directly to the remainder of the calculation
    """

    binMin = data.loc[:, prop].values.min()
    binMax = data.loc[:, prop].values.max()

    if logBin:
        if binMin*binMax < 0:
            # If the binning crosses 0, log binning will not work.
            binPos = data.loc[data.loc[:, prop] > 0, prop].values
            binNeg = data.loc[data.loc[:, prop] < 0, prop].values
            binPosMax = binPos.max()
            binPosMin = binPos.min()
            binNegMax = abs(binNeg).max()
            binNegMin = abs(binNeg).min()
            velBinPos = np.logspace(np.log10(binPosMin*0.99),
                                    np.log10(binPosMax*1.01),
                                    num=int((nBins-2)/2))
            velBinNeg = np.logspace(np.log10(binNegMin*0.99),
                                    np.log10(binNegMax*1.01),
                                    num=int((nBins-2)/2))
            velBinNeg *= -1
            velBinNeg.sort()
            velBin = np.concatenate((velBinNeg, [0], velBinPos))
        elif binMin > 0:
            velBin = np.logspace(np.log10(binMin*0.99),
                                 np.log10(binMax*1.01), num=nBins)
        elif binMin < 0:
            velBin = np.logspace(np.log10(binMax*-.99),
                                 np.log10(binMin*-1.01), num=nBins)
            velBin *= -1
    else:
        if binMin < 0:
            binMin *= 1.01
        else:
            binMin *= 0.99
        if binMax > 0:
            binMax *= 1.01
        else:
            binMax *= 0.99
        if binMin*binMax < 0:
            # If the space crosses 0, insert 0 at the interface
            velBin = np.linspace(binMin, binMax, num=nBins)
            zc = crossings_nonzero_all(velBin)+1
            velBin = np.insert(velBin, zc, [0])

        else:
            velBin = np.linspace(binMin, binMax, num=nBins)
    velVal = (velBin[1:]+velBin[:-1])/2
    normFreq, velBinCopy = np.histogram(data.loc[:, prop].values, bins=velBin,
                                        weights=data.loc[:, 'EleVol'],
                                        density=True)
    return normFreq, velVal, velBin


def extractParams(fileName, nPil=2):
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
    yInlet = data.y.values.min()
    inletData = data.loc[(data.y <= yInlet+regionSize), :]
    yOutlet = data.y.values.max()
    outletData = data.loc[(data.y >= yOutlet-regionSize)]
    inletP = inletData.p.values.mean()
    outletP = outletData.p.values.mean()
    deltaP = abs(inletP - outletP)  # Pa
    dx = abs(inletData.y.mean() - outletData.y.mean())
    return deltaP, q0, dx


def calcVortVelAngle(data, uxName, uyName, uzName, wxName, wyName, wzName):
    """each entry is just the string"""
    ux = data.loc[:, uxName].values
    uy = data.loc[:, uyName].values
    uz = data.loc[:, uzName].values
    wx = data.loc[:, wxName].values
    wy = data.loc[:, wyName].values
    wz = data.loc[:, wzName].values
    velMag = np.sqrt(ux**2+uy**2+uz**2)
    vortMag = np.sqrt(wx**2+wy**2+wz**2)
    angle = (ux*wx+uy*wy+uz*wz)/(velMag*vortMag)
    angle = np.degrees(np.arccos(angle))
    return angle


def genOutputFolderAndParams(dataDir, caseName, caseExt, nBins, logBins,
                             binProp, regionName='Result', dataRegionX=None,
                             dataRegionY=None, dataRegionZ=None):
    if logBins:
        binType = 'log'
    else:
        binType = 'linear'
    outputPath = '..\\{}-{}-{} {} bins\\'.format(regionName,
                                                 binProp, nBins, binType)
    outputFile = outputPath+'Parameters.txt'
    if not os.path.isdir(outputPath):
        os.mkdir(outputPath)
    with open(outputFile, "w") as outFile:
        outFile.write("Data directory: {}\n".format(dataDir))
        outFile.write("Case name: {}\n".format(caseName))
        outFile.write("Case extension: {}\n".format(caseExt))
        outFile.write("Binned property:{}\n".format(binProp))
        outFile.write("Number of bins: {}\n".format(nBins))
        outFile.write("Bin type: {}\n".format(binType))
        outFile.write("X Region: {}\n".format(dataRegionX))
        outFile.write("Y Region: {}\n".format(dataRegionY))
        outFile.write("Z Region: {}\n".format(dataRegionZ))
    return outputPath

# Read through files in a directory


#workingDir = "..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\RawData\\"
workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\ChemData\\RawData\\"
#workingDir = "TestData"
caseName = "TwoInletsTwoColumns_v5."
caseExt = "\.chemdata.txt$"
calcFlow = False  # Do Pressure/Flow rate fitting? Only valid with flow
writeMeta = True  # Create new metadata files
vortAng = True  # Calculate the angle between velocity and vorticity vector, will generate data column "angle"
calcChem = False  # Do calculations for PDF from chemistry

#PDF Properties

binVel = True  # True to bin velocties, false to skip
dataRegionX = [150, 350]
dataRegionY = [-550, -250]  # [-5000, 250]
regionName = 'Pillar gap'
nBins = 100
logBins = False  # True to use log spaced bins, False to use linear bins
nPil = 1  # Number of pillars in file specification
binProp = 'dCdt'  # Name of column to run PDF on, use 'angle' to do a vort./vel. angle analysis

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')

# Check for existence of a metadata file which should be something like:
# caseName+"_meta.csv"
# Purpose of metadata file is to say what we've run already
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2',
                                     'd', 'Re', 'dP', 'q', 'l'])
outFile = genOutputFolderAndParams(workingDir, caseName, caseExt,
                                   nBins, logBins, binProp,
                                   regionName=regionName,
                                   dataRegionX=dataRegionX,
                                   dataRegionY=dataRegionY)
print(outFile)
for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        # Check for fileName already in metaData, skip if so
        data = dataLoader(fileName, type=caseExt[2:-1])
        data = subSelectData(data, xRange=dataRegionX, yRange=dataRegionY)
        params = extractParams(fileName, nPil)
        params['dP'], params['q'], params['l'] = calcFlowPress(data, params)
        params['fileName'] = fileName
        if vortAng:
            data.loc[:, 'angle'] = calcVortVelAngle(data, 'u', 'v', 'w',
                                                    'vortX', 'vortY', 'vortZ')
        if calcChem:
            kVal = data.loc[data.index[0], 'k']
            dCdt = data.h2o2.values*data.tcpo.values*kVal
            data.loc['dCdt'] = dCdt
            elementVol = data.eleVol.values
            params['totalVol'] = np.sum(elementVol)
            params['totalProd'] = np.sum(np.product(data.cProduct.values, elementVol))
            params['totalTCPO'] = np.sum(np.product(data.tcpo.values, elementVol))
            params['totalH2O2'] = np.sum(np.product(data.h2o2.values, elementVol))
            params['dCdtAvg'] = np.mean(dCdt)
            params['dCdtStd'] = np.std(dCdt)
            constC = (data.tcpo.values+data.cProduct.values)  # Conservative component
            dCdtNorm = data.h2o2.values*data.tcpo.values/(1.0**2)

            params['conservative'] = data.cProduct.values*data.eleVol.values
        if binVel:
            normFreq, valMean, valBin = \
                producePDF(data, nBins=nBins, logBin=logBins, prop=binProp)
            pdfData = {'normFreq': normFreq, 'valMean': valMean,
                       'leftBin': valBin[:-1], 'rightBin': valBin[1:]}
            velDF = pd.DataFrame(pdfData)
            velDF.to_csv(outFile+fileName[:-4]+"_histogram.csv")
            plt.figure()
            plt.plot(valMean, normFreq)
            plt.title(binProp)
            plt.xlabel('Average value of bin')
            plt.ylabel('Normalized Frequency (.)')
            plt.savefig(outFile+fileName[:-4]+"_linear.png")
            plt.yscale('log')
            #plt.xscale('log')
            plt.savefig(outFile+fileName[:-4]+"_log.png")
            plt.close()
        metaData = metaData.append(params, ignore_index=True)


flowFitData = pd.DataFrame([], columns=['r1', 'r2', 'd', 'linA', 'linB',
                                        'quadA', 'quadB', 'quadC', 'expA',
                                        'expB'])
# Obtain the unique geometries used (i.e. by d, r1, and r2)
uniqueParam = metaData.loc[~metaData.duplicated(['d', 'r1', 'r2']), :]
colorPalette = sns.color_palette('deep', n_colors=len(uniqueParam.index))

f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
if calcFlow:
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
