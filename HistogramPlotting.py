import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
import seaborn as sns
from itertools import cycle

"""Purpose of script is to open and plot generated histogram files, skipping
the analysis, but instead just putting out the images.

TODO: Data that is collated should be output as a csv (or at least have a flag to output as a CSV)

CURRENTLY SWEEPS A WORKING FOLDER, TURN INTO A FUNCTION

WOULD BE NICE TO MAKE A COMPARITOR FUNCTION THAT CAN TAKE TWO FILE NAMES (OR TWO COMPARISON DIRECTORIES)

- Better better yet, have each metadata set write what folder it came from, which it can then use as a label


HEY SOME STUFF YOU NEED TO DO

-Calculate distribution median (point at which each area is divided half and half)
-create a plot that compares various aggergate statistics to a given parameter



 """

class DataSet:
    def __init__(self, workingDir, ext, label, caseName, smooth=False, window=5):
        self.workingDir = workingDir
        self.caseExt = ext
        self.label = label
        self.caseName = caseName
        self.smooth = smooth
        self.window = window
        self.dataSet = dataExtraction(self.workingDir, self.caseName,
                                      self.caseExt, self.smooth, self.window)
        self.metaData = dataSetPlot(self.dataSet)


sns.set_context('poster', font_scale=1.25)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Cambria'


def importParams(dir='.'):
    # Scans for a metafile containing the file name associated parameters
    metaPat = re.compile('.*_meta.csv')
    for fileName in os.listdir(dir):
        if re.match(metaPat, fileName):
            metaData = pd.read_csv(fileName, index_col=0)
            return metaData
        else:
            continue
    return pd.DataFrame([])


def extractParams(fileName, metaData):
    # Produces a dictionary of experimental parameters based on the file name
    # Now also will attempt to extract those params from the existing file
    # Need to update this to capture c and k values
    # Really need to update this to pull the correct file from the metaData file if available
    if not metaData.empty:
        # This will break if you are not outputting the metaData file from VelocityBinning.py
        origFileName = fileName.replace('_histogram.csv', '.txt')
        entry = metaData.loc[metaData.fileName == origFileName]
        res = entry.to_dict('records')[0]
    else:
        # Since the bining script ALWAYS outputs metadata this should not ever be called
        r1Pat = re.compile('r1_(\d+\.?\d*)_')
        r2Pat = re.compile('r2_(\d+\.?\d*)_')
        rePat = re.compile('Re(\d+\.?\d*)\.(chem|flow)data')
        dPat = re.compile('d(\d+\.?\d*)_')
        cPat = re.compile('c(\d+\.?\d*)')
        kPat = re.compile('k(\d+\.?\d*)_')
        r1Val = float(re.search(r1Pat, fileName).group(1))
        r2Val = float(re.search(r2Pat, fileName).group(1))
        dVal = float(re.search(dPat, fileName).group(1))
        reVal = float(re.search(rePat, fileName).group(1))
        cVal = float(re.search(cPat, fileName).group(1))
        kVal = float(re.search(kPat, fileName).group(1))
        vel = reVal*4.3E-7/500E-6  # Viscosity, and 500 um channel width
        daAdv = kVal*cVal/1000*2E-6*r1Val/vel  # Convert r1 to microns and give it as a diameter
        daDiff = kVal*cVal/1000*(2E-6*r1Val)**2/3E-9  # Assuming TCPO diff. coeff.
        pe = vel*r1Val*2E-6/3E-9
        res = {'r1': r1Val, 'r2': r2Val, 'd': dVal,
               'Re': reVal, 'c': cVal, 'k': kVal, 'velChar': vel, 'DaAdv': daAdv,
               'DaDiff': daDiff, 'Pe': pe}
    stokesPat = re.compile('_Stokes_')
    if re.search(stokesPat, fileName):
        res['Flow'] = 'Stokes'
    else:
        res['Flow'] = 'NS'
    rePillar = res['Re']*res['r1']/250  # Assume 500 um channel width
    res['RePil'] = rePillar
    return res


def dataExtraction(workingDir, caseName, caseExt, smooth=False, window=5):
    os.chdir(workingDir)
    filePat = re.compile(caseName+'.*?'+caseExt)
    fileList = os.listdir('.')
    dataSets = {}
    for fileName in fileList:
        if re.match(filePat, fileName):
            # params = extractParams(fileName)
            data = pd.read_csv(fileName, header=0)
            if smooth:
                # Smooth data with rolling average of size window
                data = dataSmoothing(data, window)
                # Drop NaNs
                data = data.dropna()
            dataSets[fileName] = data
    return dataSets


def dataSetPlot(dataSets, metaData=pd.DataFrame([]), linestyle='-', smooth=0):
    # In this case, key refers to the filename for the given data set
    workingDirParams = importParams()
    # print(workingDirParams
    if not metaData.empty:
        metaData = pd.DataFrame([], columns=['r1', 'r2', 'd', 'Re', 'Flow', 'PDFmean', 'PDFstd'])
    for key in dataSets:
        print(key)
        data = dataSets[key]
        params = extractParams(key, workingDirParams)
        dataMean, dataVar = pdfStats(data)
        params['PDFmean'] = dataMean
        params['PDFstd'] = np.sqrt(dataVar)
        metaData = metaData.append(params, ignore_index=True)
        # dataMid = np.mean([np.min(data.normFreq), np.max(data.normFreq)])
        # pmf, xPmf = genPMF(data)
        a1 = ax1.plot(data.valMean, data.normFreq,
                      label=key+'smooth {}'.format(smooth), ls=linestyle)
        c = a1[0].get_color()
        # ax1.plot([dataMean, dataMean],
        #          [np.min(data.normFreq), np.max(data.normFreq)],
        #          ls='-', color='k')
        # ax1.plot([dataMean-np.sqrt(dataVar), dataMean+np.sqrt(dataVar)],
        #          [dataMid, dataMid], ls='--', color='k')

        a2 = ax2.plot(data.valMean,
                      data.normFreq, label=key+'smooth {}'.format(smooth),
                      ls=linestyle)
        c = a2[0].get_color()
        # ax2.plot([dataMean, dataMean],
        #          [np.min(data.normFreq), np.max(data.normFreq)],
        #          ls='--', color='k')
        # ax2.plot([dataMean-np.sqrt(dataVar), dataMean+np.sqrt(dataVar)],
        #          [dataMid, dataMid], ls='--', color='k')
        diffPDF = np.diff(data.normFreq)/np.diff(data.valMean)
        ax3.plot(data.valMean[1:], diffPDF, label=key, ls=linestyle)
        # ax5.plot(xPmf, pmf, label=key, ls=linestyle)
# HEY UNIFY THE NAMING IN THE MAIN SCRIPT SINCE IT'S NOT ALWAYS VELOCITIES
    return metaData


def dataSmoothing(data, window=5):
    dataSmooth = data.rolling(window, center=True).mean()
    return dataSmooth


def pdfStats(data):
    dx = data.rightBin-data.leftBin
    eVal = np.sum(data.valMean*data.normFreq*dx)
    varVal = np.sum(data.normFreq*(data.valMean-eVal)**2*dx)
    return eVal, varVal


def genGaussian(mu, variance):
    sigma = np.sqrt(variance)
    x = np.linspace(mu-5*sigma, mu+5*sigma, 1000)
    y = stats.norm.pdf(x, mu, sigma)
    return x, y


def semilogYFitting(data, xPropName, yPropName, xRange):
    subData = data.loc[(data.loc[:, xPropName] >= xRange[0]) & (data.loc[:, xPropName] <= xRange[1]), :]
    xData = subData.loc[:, xPropName].values
    yData = np.log10(subData.loc[:, yPropName].values)
    fit = np.polyfit(xData, yData, 1)
    yEst = np.power(10, np.polyval(fit, xData))
    return fit, xData, yEst


def metaPlot(metaData, prop='Re', flowCond='NS', label=None, marker='o'):
    if not label:
        label = flowCond+" "+prop
    subData = metaData.loc[metaData.loc[:, 'Flow'] == flowCond, :]
    subData.sort_values(by=prop, inplace=True)
    ax4.errorbar(subData.loc[:, prop], subData.loc[:, 'PDFmean'],
                 yerr=subData.loc[:, 'PDFstd'], ls='-', marker=marker,
                 capsize=2, label=label)
    ax4.set_xlabel(prop)
    ax6.plot(subData.loc[:, 'RePil'], subData.loc[:, prop],
             ls='-', marker=marker,
             label=label)
    ax6.set_ylabel(prop)
    ax6.set_xlabel('Re')
    maxVal = subData.loc[:, 'PDFmean'].max()
    ax7.plot(subData.loc[:, prop], subData.loc[:, 'PDFmean']/maxVal,
                 ls='-', marker=marker,
                 label=label+' Max val: {}'.format(maxVal))
    ax7.set_xlabel(prop)
    ax7.set_ylabel('Max normalized')
    ax7.set_title('Watch out for what the max means')
    ax7.set_ylim([0.3, 1.05])
    ax7.set_xlim([-1, 105])


    maxVal = subData.loc[:, 'volWeightedMean'].max()
    ax8.plot(subData.loc[:, prop], subData.loc[:, 'volWeightedMean']/maxVal,
                 ls='-', marker=marker,
                 label=label+' Max val: {}'.format(maxVal))
    ax8.set_xlabel(prop)
    ax8.set_ylabel('Volume Weighted Mean Normalized to Max (.)')
    ax8.set_title('Watch out for what the max means')
    ax8.legend(loc=0)
    ax8.set_ylim([0.3, 1.05])
    ax8.set_xlim([-1, 105])


    ax9.plot(subData.loc[:, prop], subData.loc[:, 'volWeightedMean'],
             ls='-', marker=marker, label=label)
    ax9.set_xlabel(prop)
    ax9.set_ylabel('Volume Weighted Mean (.)')
    ax9.legend(loc=0)
    ax9.set_xlim([-1, 105])
    return


def genPMF(data):
    """Perform trapz calculation manually so we can use cumsum() to get the PMF
    Not saving it in the main dataframe since we have to reduce the dimensionality
    to estimate the integral

    PMF is technically defined for a given region (i.e. from where to where do
    I integrate the PDF?)

    For this calculation, I produce a result, g(x) = integral from 0 to x of
    f(a)da using the trapezoidal rule
    """
    f = data.normFreq.values
    x = data.valMean.values
    fInt = (f[1:]+f[:-1])/2*(x[1:]-x[:-1])
    PMF = np.cumsum(fInt)
    xVal = (x[1:]+x[:-1])/2  # Use average to give what x is the corresponding value
    return PMF, xVal


def runDataSet(workingDir, caseName, caseExt, label, smooth, window, metaData, props):
    """Helper function which loads data, produces the associated metaData, then plots it
    Would be nice to maybe pack these things up (workingDir, caseName, ext, Label)
    all belong together"""
    dataSet = dataExtraction(workingDir, caseName, caseExt, smooth, window)
    metaData = dataSetPlot(dataSet, metaData, smooth=window, linestyle='-')
    for prop in props:
        metaPlot(metaData, prop=prop, label=label+' '+prop)
    return dataSet, metaData


quickVol = {0.1: 1.48E-13, 1: 6.02E-14, 10: 5.85E-14, 50: 4.58E-13}
#TODO: Change this move to the base directory containing the data, then operate out of that folder
window = 5
smooth = False
fitRange = np.array([85, 90])
prop = 'RePil'  # Options include: Re, d, RePil, DaAdv, DaDiff, Pe, reactorConserv
prop2 = None # #'posMRT' # Lets you plot multiple properties vs Re, beware axis scaling
# fitRange = np.array([65, 85])
# workingDirA = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\FlowData\\RecircZoneBasic-velMag-100 linear bins\\"
#workingDirA = "..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\Pillar Gap-angle-180 linear bins"
workingDirA = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\Bottom gap pillar exclusive-constC-100 linear bins\\"
#workingDirA = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\Pillar gap pillar exclusive-constC-100 linear bins\\"
# workingDir = "."
caseNameA = "TwoPillar_v6_ExF_"
caseNameA = "TwoPillar_v6_ExF_c3_k2000_"
caseExtA = "d25_Re.*\.chemdata_histogram\.csv"
# caseExtA = "_Re.*\.chemdata_histogram\.csv"
# labelA = "Pillar Inclusive"
labelA = "25 um-tracer"
# workingDirB = "..\\..\\..\\..\\..\\Multipillar\\Normal\\FlowData_Normal\\200 log bins - 250 to -2500"
#workingDirB = "..\\..\\..\\..\\..\\..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\Pillar gap-angle-180 linear bins"
workingDirB = "..\\Bottom gap pillar exclusive-tcpo-100 linear bins\\"
#workingDirB = "."
caseNameB = "TwoPillar_v6_ExF_c3_k2000_"
caseExtB = "d25_Re.*\.chemdata_histogram\.csv"
labelB = "25 um-TCPO"

# workingDirC = "..\\Bottom gap pillar exclusive-h2o2-100 linear bins\\"
# caseNameC = "TwoPillar_v6_ExF_c3_k2000_"
# caseExtC = "d100_Re.*\.chemdata_histogram\.csv"
# labelC = "H2O2"


# Plot for everything
f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f3, ax3 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f4, ax4 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f5, ax5 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f6, ax6 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f7, ax7 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f8, ax8 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f9, ax9 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))


metaData = pd.DataFrame([], columns=['r1', 'r2', 'd', 'Re', 'Flow', 'PDFmean', 'PDFstd'])
dataSetA = dataExtraction(workingDirA, caseNameA, caseExtA, smooth, window)
metaDataA = dataSetPlot(dataSetA, metaData, smooth=window, linestyle='-')

dataSetB = dataExtraction(workingDirB, caseNameB, caseExtB, smooth, window)
metaDataB = dataSetPlot(dataSetB, metaData, smooth=window,linestyle='-')

# dataSetC = dataExtraction(workingDirC, caseNameC, caseExtC, smooth, window)
# metaDataC = dataSetPlot(dataSetC, metaData, smooth=window,linestyle='-')

markerCycle = cycle(['o', '^', 's', 'd', 'D'])
metaPlot(metaDataA, prop=prop, flowCond='NS', label=labelA+' '+prop, marker=next(markerCycle))
metaPlot(metaDataB, prop=prop, flowCond='NS', label=labelB+' '+prop, marker=next(markerCycle))
# metaPlot(metaDataC, prop=prop, flowCond='NS', label=labelC+' '+prop, marker=next(markerCycle))

if prop2:
    metaPlot(metaDataA, prop=prop2, flowCond='NS', label=labelA+' '+prop2)
    # metaPlot(metaDataB, prop=prop2, flowCond='NS', label=labelB+' '+prop2)
    # metaPlot(metaDataC, prop=prop2, flowCond='NS', label=labelC+' '+prop2)
ax1.set_title("PDFs")
ax2.set_title("PDFs")
ax3.set_title("Differentiated PDF")
ax1.set_xlabel("Value")
ax1.set_ylabel("Normalized freq.")
ax2.set_xlabel("Value")
ax2.set_ylabel("Normalized freq.")
# ax1.legend(loc=0) # Legends get too long for this
ax1.set_yscale('log')
# ax2.legend(loc=0)
# ax3.legend(loc=0)

ax4.set_ylabel('Mean of PDF')
ax3.set_yscale('log')
# plt.xscale('log')
ax4.legend(loc=0)
ax5.set_title("PMFs")
ax5.set_xlabel('Value')
ax5.set_ylabel('PMF')
ax5.legend(loc=0)
ax5.set_yscale('log')
ax6.legend(loc=0)

sns.despine(f1)
sns.despine(f2)
sns.despine(f3)
sns.despine(f4)
sns.despine(f5)
sns.despine(f6)
sns.despine(f7)
sns.despine(f8)
sns.despine(f9)

plt.ion()
plt.show()
