import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns


"""TODO: Convert script inputs into user prompts, see input(). Will need to set up to take default values as well."""


def stringVal(pattern, string):
    """Given a regular expression searching for a single instance of a numeric
    value, function will return a float of that value, or -1 if no match is found."""
    val = re.search(pattern, string)
    if not val:
        return -1
    else:
        return float(val.group(1))


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


def extractParams(fileName, nPil=2):
    # Produces a dictionary of experimental parameters
    rePat = re.compile('Re(\d+\.?\d*).')
    dPat = re.compile('d(\d+\.?\d*)_')
    cPat = re.compile('c(\d+\.?\d*)')
    kPat = re.compile('k(\d+\.?\d*)_')
    dVal = stringVal(dPat, fileName)
    reVal = stringVal(rePat, fileName)
    cVal = stringVal(cPat, fileName)
    kVal = stringVal(kPat, fileName)
    if nPil == 2:
        r1Pat = re.compile('r1_(\d+?)_')
        r2Pat = re.compile('r2_(\d+?)_')
        r1Val = stringVal(r1Pat, fileName)
        r2Val = stringVal(r2Pat, fileName)
        res = {'r1': r1Val, 'r2': r2Val, 'd': dVal,
               'Re': reVal, 'c': cVal, 'k': kVal}
    if nPil == 1:
        rPat = re.compile('r(\d+?)_')
        rVal = stringVal(rPat, fileName)
        res = {'r1': rVal, 'r2': rVal, 'd': dVal,
               'Re': reVal, 'c': cVal, 'k': kVal}
    return res


def pillarGapCalculation(r1, r2, d, includePillar=True):
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
    if includePillar:
        y1 = yPil+d/2+r1 # Includes the pillar itself, allowing the box to cover volume near the pillar far away from the centerline
        y2 = yPil-d/2-r2
    else:
        y1 = yPil+d/2
        y2 = yPil-d/2
    return [x1, x2], [y1, y2]


def loadData(fileName):
    """Read in data from a single fileName, determine the relevant
    model parameters from the filename, strip out any data points outside the
    "recirculation" zone, then return the stripped data and the parameters.

    Is there a reason we would prefer to import the data a different way?
    Perhaps reading in the data directly using a file reader?
    """
    data = pd.read_table(fileName, sep='\s+', skiprows=8,
                         names=['x', 'y', 'z', 'sID', 'velMag'])
    params = extractParams(fileName, nPil=1)
    params['fileName'] = fileName
    xRange, yRange = pillarGapCalculation(params['r1'], params['r2'], params['d'])
    data = subSelectData(data, xRange=xRange, yRange=yRange)
    return data, params


def groupStreamlines(data):
    xVal = data.x.values
    yVal = data.y.values
    zVal = data.z.values
    sID = data.sID.values
    dist = np.sqrt(np.square(xVal[:-1]-xVal[1:])+np.square(yVal[:-1]-yVal[1:])+np.square(zVal[:-1]-zVal[1:]))*1E-6
    sIDMask = sID[1:] == sID[:-1]  # Find the points where the switch occurs
    dist = dist*sIDMask  # Sets the values where the streamlines change to 0
    # Make a mask which sets dist to 0 where the next element is on a new streamline
    dist = np.append(dist, 0)  # The last point will contribute no travel time, also aligning the data
    data['dist'] = dist  # Attach distance to a given coordinate (note distance is in m)
    # ID indices where the streamline number changes, set dist to 0 at those points
    data['travelTime'] = dist/data.velMag  # Calculate travel time
    sIDGroups = data.groupby(data.sID)  # Group everything by streamline ID
    return sIDGroups


def calcPDFStats(bins, freq):
    """ Calcuates the mean and variance of a histogram produced by
    np.histogram().

    inputs:

    bins: N bin edges
    freq: N-1 frequency values corresponding to each bin

    returns:

    mean: the mean of the PDF
    variance: the pdf variance
    """
    rightBin = bins[1:]
    leftBin = bins[:-1]
    binMean = (leftBin+rightBin)/2
    dx = rightBin
    mean = np.sum(binMean*freq*dx)
    variance = np.sum(freq*(binMean-mean)**2*dx)
    return mean, variance


def generatePDFs(workingDir, caseName, ext, testMode, nBins, outputPath):
    """Generates PDFs and outputs both the histograms, histogram plots,
    and statistical information about the PDFs.

    For now, log(time) and streamline length are binned.

    INPUTS
    workingDir: Working directory where data to be binned is stored
    caseName: Regex expression to produce PDFs using only certain files
    ext: Extension for the PDFs of interest
    testMode: Flag which will run only one file and then return.

    """
    # Intialize output folder and metaData data structure

    filePat = re.compile('(.*)('+ext+')')
    fileList = os.listdir('.')
    metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2', 'd', 'Re',
                                         'c', 'k'])
    for f in fileList:
        if re.match(filePat, f):
            print(f)
            # Load in data set and group by streamline
            dataSet, params = loadData(f)
            gStream = groupStreamlines(dataSet)
            # Calculate travel times
            travelTime = gStream.travelTime.sum().values
            streamLen = gStream.dist.sum().values
            if logVal:
                travelTime = np.log10(travelTime)
                streamLen = np.log10(streamLen)
            # refTime = streamLen/gStream.velMag.mean().values
            timePDF, timeBins = np.histogram(travelTime, bins=nBins,
                                             density=True)

            # timeBins = np.power(10, timeBins)  # Convert back to non logspace times
            lenPDF, lenBins = np.histogram(streamLen, bins=nBins,
                                           density=True)

            output = pd.DataFrame({'lenLBins': lenBins[:-1],
                                   'lenRBins': lenBins[1:],
                                   'lenBinMean': (lenBins[:-1]+lenBins[1:])/2,
                                   'lenPDF': lenPDF,
                                   'timeLBins': timeBins[:-1],
                                   'timeRBins': timeBins[1:],
                                   'timeBinMean': (timeBins[:-1]+timeBins[1:])/2,
                                   'timePDF': timePDF})
            params['lenPDFMean'], params['lenPDFVar'] = calcPDFStats(lenBins, lenPDF)
            params['timePDFMean'], params['timePDFVar'] = calcPDFStats(timeBins, timePDF)
            params['timeMean'] = np.mean(travelTime)
            params['timeStd'] = np.std(travelTime)
            params['lenMean'] = np.mean(streamLen)
            params['lenStd'] = np.std(streamLen)
            params['RePil'] = params['Re']*params['r1']*2/500  # RePil = r1*2/500*Re
            metaData = metaData.append(params, ignore_index=True)
            output.to_csv(outputPath+f+".histogram.csv")
            plt.figure(1)
            plt.plot((timeBins[1:]+timeBins[:-1])/2, timePDF)
            plt.title('timePDF')
            if logVal:
                plt.xlabel('log(Time)')
            else:
                plt.xlabel('Time (s)')
            plt.ylabel('PDF')
            plt.yscale('log')
            plt.savefig(outputPath+f[:-4]+'_timePDF.png', dpi=600)
            sns.despine()
            plt.figure(2)
            plt.plot((lenBins[1:]+lenBins[:-1])/2, lenPDF)
            plt.title('lengthPDF')
            if logVal:
                plt.xlabel('log(Streamline Length)')
            else:
                plt.xlabel('Streamline Length (m)')
            plt.ylabel('PDF')
            sns.despine()
            plt.savefig(outputPath+f[:-4]+'_lenPDF.png', dpi=600)
            if testMode:
                plt.ion()
                fig3 = plt.figure(3)
                ax3 = fig3.add_subplot(111, projection='3d')
                for sID in gStream.sID.min().values:
                    streamLine = dataSet.loc[dataSet.sID == sID, :]
                    ax3.plot(streamLine.x, streamLine.y, streamLine.z)
                plt.title('Streamlines')
                return metaData
            plt.close('all')
    metaData.to_csv(outputPath+caseName+'_metaData.csv')
    return metaData


def generateMeanPlots(metaData, logVal):
    # Should correctly consider if log of values or not.
    f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
    f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))

    f3, ax3 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
    f4, ax4 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))

    for d in metaData.d.unique():
        subData = metaData.loc[metaData.d == d, :]
        if logVal:
            x = subData.RePil
            yTime = np.power(10, subData.timeMean)
            yerrTime = np.power(10, subData.timeStd)
            yLen = np.power(10, subData.lenMean)
            yerrLen = np.power(10, subData.lenStd)
        else:
            x = subData.RePil
            yTime = subData.timeMean
            yerrTime = np.sqrt(subData.timeStd)
            yLen = subData.lenMean
            yerrLen = np.sqrt(subData.lenStd)

        ax1.errorbar(x, yTime, yerr=yerrTime,
                     ls='None', marker='.', label=float(d))
        ax2.errorbar(x, yLen, yerr=yerrLen,
                     ls='None', marker='.', label=float(d))
        ax3.plot(x, yTime, ls='None', marker='.', label=float(d))
        ax4.plot(x, yLen, ls='None', marker='.', label=float(d))

    ax1.set_xlabel('RePil')
    ax2.set_xlabel('RePil')
    ax3.set_xlabel('RePil')
    ax4.set_xlabel('RePil')
    ax1.set_ylabel('Mean Residence Time (s)')
    ax2.set_ylabel('Mean Streamline Length (m)')
    ax3.set_ylabel('Mean Residence Time (s)')
    ax4.set_ylabel('Mean Streamline Length (m)')
    ax1.set_title('Streamline Residence Time')
    ax2.set_title('Streamline Length')
    ax3.set_title('Streamline Residence Time')
    ax4.set_title('Streamline Length')
    sns.despine(f1)
    sns.despine(f2)
    sns.despine(f3)
    sns.despine(f4)
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax4.legend()
    return ax1, ax2, ax3, ax4


def plotPDF():
    f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
    f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
    return ax1, ax2
"""SCRIPT INPUTS"""


workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\RecircZoneStreamlines\\PillarGap\\"
caseName = "TwoPillar_"
ext = ".velStreamline.txt"
testMode = False  # Runs on one file, produces plots, then stops PDF calculation
nBins = 100  # Number of bins to use for PDF
calculatePDFs = False  # Flag to toggle calculation of PDFs.
logVal = False  # Bin log values instead of the actual values

"""MAIN SCRIPT"""

sns.set_context('talk')
plt.rcParams['svg.fonttype'] = 'none'
os.chdir(workingDir)
if logVal:
    outputPath = ".\\log value PDF {} bins\\".format(nBins)
else:
    outputPath = ".\\PDF {} bins\\".format(nBins)
if not os.path.isdir(outputPath):
    os.mkdir(outputPath)
if calculatePDFs:
    metaData = generatePDFs(workingDir, caseName, ext, testMode, nBins,
                            outputPath)
else:
    metaData = pd.read_csv(outputPath+caseName+"_metaData.csv")

ax1, ax2, ax3, ax4 = generateMeanPlots(metaData, logVal)
plt.ion()
plt.show()
