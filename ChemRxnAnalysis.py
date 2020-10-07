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


def producePDF(data, dataCol, nBins=1000, logBin=True):
    """
    INPUT:
    data: data output from comsol, must include "eleVol" col and column
    named in dataCol
    dataCol: String of the named column to bin
    nBins: Number of bins you want to use, default is 1000
    logBin: True-> bins are evenly spaced in log space,
    otherwise, linear spacing

    OUTPUT
    normFreq: Normalized frequencies for the defined bins by bin size
    meanVal: The mean velocity in each bin
    groups: The binned data group
    bins: Velocity bin limits used
    Produce a velocity PDF from the input data

    Option, bin velocities such that the mean and median are within
    some tolerance or the RSD meets a certain criteria

    Max velocity is taken as 1.01*max due to the way np.digitize works.
    For a linear bin, without doing this, the last value will not belong to a
    bin that has a proper size definition.
    """
    if logBin:
        bins = np.logspace(np.log10(data.loc[:, dataCol].min()),
                           np.log10(data.loc[:, dataCol].max()*1.01),
                           num=nBins)
    else:
        bins = np.linspace(data.loc[:, dataCol].min(),
                           data.loc[:, dataCol].max()*1.01,
                           num=nBins)
    binsSize = bins[1:]-bins[:-1]
    data.loc[:, 'binID'] = np.digitize(data.loc[:, dataCol], bins)
    groups = data.groupby(data.binID)
    meanVal = groups[dataCol].mean()
    # Weight frequencies by included volume and normalize to bin size
    weightedFreq = groups.eleVol.sum()*groups.size() \
        / groups.binID.median().apply(lambda x: binsSize[x-1])
    # Calculate area under freq-vel curve to obtain a PDF
    totalArea = np.trapz(weightedFreq, x=meanVal)
    normFreq = weightedFreq/totalArea
    return normFreq, meanVal, groups, bins

# Read through files in a directory


workingDir = "..\\Comsol5.4\\Multipillar\\Normal\\ChemRxn_Normal\\First pillar gap\\"
# workingDir = "."
caseName = "Multipillar_v5.2_Normal_r50_d50"
caseExt = "\.chemdata\.txt$"
writeMeta = True
binProd = True
binRate = True
dataRegionY = [-150, -100]
dataRegionX = [200, 300]
nBins = 100

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
        data = subSelectData(data, xRange=dataRegionX, yRange=dataRegionY)
        data['cProductNorm'] = data.cProduct/1.0
        if binProd:
            prodFreq, prodMean, prodGroups, prodBins = producePDF(data, 'cProductNorm', nBins=nBins, logBin=False)
            prodData = {'normFreq': prodFreq, 'valMean': prodMean}
            velPDF = pd.DataFrame(prodData)
            velPDF.to_csv(fileName[:-4]+"_ProductHistogram.csv")
            plt.figure()
            plt.plot(prodMean, prodFreq)
            plt.xlabel('Average value')
            plt.ylabel('Normalized Frequency (.)')
            plt.savefig(fileName[:-4]+"_Product_linear.png")
            plt.yscale('log')
            # plt.xscale('log')
            plt.savefig(fileName[:-4]+"_Product_log.png")
            plt.close()
        data['dCdt'] = data.h2o2*data.tcpo*data.k
        data['dCdtNorm'] = data.h2o2*data.tcpo/(1.0**2)
        if binRate:
            rateFreq, rateMean, rateGroups, rateBins = producePDF(data, 'dCdtNorm', nBins=nBins, logBin=False)
            rateData = {'normFreq': rateFreq, 'valMean': rateMean}
            velPDF = pd.DataFrame(rateData)
            velPDF.to_csv(fileName[:-4]+"_RateHistogram.csv")
            plt.figure()
            plt.plot(rateMean, rateFreq)
            plt.xlabel('Average value')
            plt.ylabel('Normalized Frequency (.)')
            plt.savefig(fileName[:-4]+"_Rate_linear.png")
            plt.yscale('log')
            # plt.xscale('log')
            plt.savefig(fileName[:-4]+"_Rate_log.png")
            plt.close()

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
        res['M'] = data.M.sum()  # Is supposed to represent quantity used to calculate the scalar dissipiation rate

        # data = subSelectData(data, yRange=dataRegionY)
        metaData = metaData.append(res, ignore_index=True)

if writeMeta:
    metaData.to_csv(caseName+"_meta.csv")
