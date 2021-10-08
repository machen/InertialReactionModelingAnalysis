import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns


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
    dist = np.sqrt(np.square(xVal[:-1]-xVal[1:])+np.square(yVal[:-1]-yVal[1:])+np.square(zVal[:-1]-zVal[1:]))*1E-6
    dist = np.append(dist, 0) # The last point will contribute no travel time
    data['dist'] = dist  # Attach distance to a given coordinate (note distance is in m)
    data['travelTime'] = dist/data.velMag  # Calculate travel time
    sIDGroups = data.groupby(data.sID)  # Group everything by streamline ID
    return sIDGroups


"""SCRIPT INPUTS"""

workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\RecircZoneStreamlines"
caseName = "TwoPillar_"
ext = ".velStreamline.txt"
testMode = True
nBins = 50

os.chdir(workingDir)
filePat = re.compile('(.*)('+ext+')')
fileList = os.listdir('.')

for f in fileList:
    if re.match(filePat, f):
        print(f)
        # Load in data set and group by streamline
        dataSet, params = loadData(f)
        gStream = groupStreamlines(dataSet)
        # Calculate travel times
        travelTime = gStream.travelTime.sum().values
        streamLen = gStream.dist.sum().values
        refTime = streamLen/gStream.velMag.mean().values
        timePDF, timeBins = np.histogram(np.log(travelTime), bins=nBins, density=True)
        lenPDF, lenBins = np.histogram(streamLen, bins=nBins, density=True)
        if testMode:
            plt.ion()
            plt.figure(1)
            plt.plot((timeBins[1:]+timeBins[:-1])/2,timePDF)
            plt.title()
            plt.figure(2)
            plt.plot((lenBins[1:]+lenBins[:-1])/2,lenPDF)
            fig3 = plt.figure(3)
            ax3 = fig3.add_subplot(111, projection='3d')
            for sID in gStream.sID.min().values:
                streamLine = dataSet.loc[dataSet.sID==sID,:]
                ax3.plot(streamLine.x,streamLine.y,streamLine.z)
            break





