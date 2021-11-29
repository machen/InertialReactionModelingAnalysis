import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
import seaborn as sns

""" TODO: THIS NEEDS TO LOAD THE METADATA AS WELL
THAT WAS GENERATED FOR EACH CASE

"""


class DataSet:
    def __init__(self, workingDir, caseName, caseExt, d,
                 label, smooth=False, window=5):
        self.workingDir = workingDir
        self.caseName = caseName
        self.caseExt = caseExt
        self.pillarGap = d
        self.label = label
        self.dataSet, self.metaData = dataExtraction(workingDir,
                                                     caseName, caseExt,
                                                     smooth, window)

    def getDataSet(self):
        return self.dataSet

    def getWorkingDir(self):
        return self.workingDir

    def getCaseName(self):
        return self.caseName

    def getCaseExt(self):
        return self.caseExt

    def getPillarGap(self):
        return self.pillarGap

    def getLabel(self):
        return self.label

    def getMetaData(self):
        return self.metaData

    def plotDataSet(self, linestyle='-'):
        for key in self.dataSet:
            data = self.dataSet[key]
            # Key is going to be a .csv file contianing the PDF data
            # We need to match it to the metadata associated
            extRe = re.compile('()('+self.caseExt[2:]+')')
            keyBase = re.search(extRe, key).group(1)+'.tif'




def flowRateConversion(q, width, height, charLen, nu=0.45):
    # Base nu is in mm2/s, so you should report this in mm and seconds as units
    reyn = q/width/height*charLen/nu
    return reyn


def extractParams(fileName):
    # Produces a dictionary of experimental parameters
    qPat = re.compile('(\d*\.?\d*)q')
    repPat = re.compile('(\d{0,3}).nd2')
    qVal = re.search(qPat, fileName).group(1)
    repVal = re.search(repPat, fileName).group(1)
    if not repVal:
        repVal = 0
    res = {'q': float(qVal), 'replicate': int(repVal)}
    return res


def loadMetaData(metaFile):
    params = pd.read_csv(metaFile, index_col=1)
    return params


def dataExtraction(workingDir, caseName, caseExt, smooth=False, window=5):
    os.chdir(workingDir)
    filePat = re.compile(caseName+'.*?'+caseExt)
    fileList = os.listdir('.')
    dataSets = {}
    for fileName in fileList:
        if re.match(filePat, fileName):
            print(fileName)
            # params = extractParams(fileName)
            data = pd.read_csv(fileName, header=0)
            if smooth:
                # Smooth data with rolling average of size window
                data = dataSmoothing(data, window)
                # Drop NaNs
                data = data.dropna()
            dataSets[fileName] = data
    metaData = loadMetaData('_meta.csv')
    return dataSets, metaData


def dataSetPlot(dataSets, metaData, d, linestyle='-', smooth=0):
    for key in dataSets:
        data = dataSets[key]
        dataMean, dataVar = pdfStats(data)
        params = extractParams(key)
        params['d'] = d
        params['fileName'] = key
        params['PDFmean'] = dataMean
        params['PDFstd'] = np.sqrt(dataVar)
        # metaData.loc[key, 'd'] = d
        # metaData.loc[key, 'PDFmean'] = dataMean
        # metaData.loc[key, 'PDFstd'] = np.sqrt(dataVar)
        metaData = metaData.append(params, ignore_index=True)
        a1 = ax1.plot(data.valMean, data.normFreq,
                      label=key+'smooth {}'.format(smooth), ls=linestyle)
        c = a1[0].get_color()
        # ax1.plot([dataMean, dataMean],
        #          [np.min(data.normFreq), np.max(data.normFreq)],
        #          ls='-', color='k')
        # ax1.plot([dataMean-np.sqrt(dataVar), dataMean+np.sqrt(dataVar)],
        #          [dataMid, dataMid], ls='--', color=c)

        a2 = ax2.plot(data.valMean,
                      data.normFreq, label=key+'smooth {}'.format(smooth),
                      ls=linestyle)
        c = a2[0].get_color()
        # ax2.plot([dataMean, dataMean],
        #          [np.min(data.normFreq), np.max(data.normFreq)],
        #          ls='-', color='k')
        # ax2.plot([dataMean-np.sqrt(dataVar), dataMean+np.sqrt(dataVar)],
        #          [dataMid, dataMid], ls='--', color=c)
        diffPDF = np.diff(data.normFreq)/np.diff(data.valMean)
        ax3.plot(data.valMean[1:], diffPDF, label=key, ls=linestyle)
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


def metaPlot(metaData, prop='q'):
    # Change me to work with the image histograms and metadata
    if prop == 'Re':
        metaData.loc[:, prop] = flowRateConversion(metaData.q/60.0, 500E-3, 100E-3, 500E-3)
    elif prop == 'ReP':
        metaData.loc[:, prop] = flowRateConversion(metaData.q/60.0, 500E-3, 100E-3, 200E-3)
    metaData.sort_values(by=prop)
    traceList = metaData.d.unique()
    for val in traceList:
        subData = metaData.loc[metaData.d == val, :]
        meanData = subData.groupby(prop).mean()
        stdData = subData.groupby(prop).std()
        ax4.errorbar(subData.loc[:, prop], subData.loc[:, 'PDFmean'],
                     yerr=subData.loc[:, 'PDFstd'], ls='none', marker='o',
                     capsize=2, label=val)
        ax5.plot(subData.loc[:, prop], subData.loc[:, 'PDFmean'],
                 ls='none', marker='o', label=val)
        maxLoc = meanData.loc[meanData.PDFmean==meanData.PDFmean.max(), 'q'].values
        ax6.errorbar(meanData.index, meanData.loc[:,'PDFmean'],
                     yerr=stdData.loc[:,'PDFstd'], ls='none', marker='o',
                     capsize=2, label= "{} max val q: {}".format(val,float(maxLoc)))
        maxVal = meanData.PDFmean.max()
        ax7.errorbar(meanData.index, meanData.loc[:,'PDFmean']/maxVal,
                     yerr=stdData.loc[:, 'PDFstd']/maxVal, ls='none', marker='o',
                     capsize=2, label= "{} max val: {}".format(val, maxVal))
        # maxVal = meanData.meanInt.max()
        # ax8.errorbar(meanData.index, meanData.meanInt/maxVal,
        #              yerr=meanData.stdInt/maxVal, ls='none',
        #              marker='o', capsize=2, label = val)
    ax4.set_xlabel(prop)
    ax4.set_ylabel('Mean of PDF')
    ax4.legend(loc=0)
    ax5.set_xlabel(prop)
    ax5.set_ylabel('Mean of PDF')
    ax5.legend(loc=0)
    ax6.set_xlabel(prop)
    ax6.set_ylabel('Mean of PDF')
    ax6.set_title('Averaged Results')
    ax6.legend(loc=0)
    ax7.set_xlabel(prop)
    ax7.set_ylabel('Mean of PDF')
    ax7.set_title('Max normalized - check max val')
    ax7.legend(loc=0)
    ax7.set_ylim([-0.1, 1.1])
    ax7.set_xlim([-1, 105])
    ax8.legend(loc=0)
    ax8.set_xlabel(prop)
    ax8.set_ylabel('Mean Intensity Normalized to Max Observed')
    ax8.set_title('Max normalized - check meaning of max val')
    return

sns.set_context('talk')
plt.rcParams['svg.fonttype'] = 'none'
smooth = False
window = 10

# Might be nice to do some averaging of lines that have the same experiemntal condition

workingDirA = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\2021-10-20-Chemilum-100um\\A3-100um\\Raw Image Pillar Gap 50 bins\\"
# workingDirB = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\2021-10-20-Chemilum-100um\\A2-100um\\Raw Image Pillar Gap 50 bins\\"
# workingDirC = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\2021-10-05-Chemilum-100um\\100 um Pillar Gap\\Raw Image Pillar Gap 50 bins\\"
# workingDirD = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\2021-11-19-Chemilum-25um\\2PD3_A4\\Raw Image Pillar Gap 50 bins\\"
workingDirE = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\2021-11-18-Chemilum-25um\\2PD3_A2\\Raw Image Pillar Gap 50 bins\\"
os.chdir(workingDirA)
caseNameA = ''
caseExtA = ".*_dark_hist\.csv"

# You must set these to the correct pillar gaps of the experiment
dA = "2021-10-20 Device A3-100um"
dB = "2021-10-20 Device A2-100um"
dC = "2021-10-05 100 um"
dD = "2021-11-19 25 um"
dE = "2021-11-18 25 um"

f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f3, ax3 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f4, ax4 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f5, ax5 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f6, ax6 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f7, ax7 = plt.subplots(1, 1, sharex='col', figsize=(12, 10)) # Plot for normalizing to max average value
f8, ax8 = plt.subplots(1, 1, sharex='col', figsize=(12, 10)) # Plot using meanIntensity rather than mean of the PDF

metaData = pd.DataFrame([], columns=['q', 'replicate', 'PDFmean', 'PDFstd'])
dataSetA, metaDataA = dataExtraction(workingDirA, caseNameA, caseExtA, smooth, window)
# dataSetB, metaDataB = dataExtraction(workingDirB, caseNameA, caseExtA, smooth, window)
# dataSetC, metaDataC = dataExtraction(workingDirC, caseNameA, caseExtA, smooth, window)
# dataSetD, metaDataD = dataExtraction(workingDirD, caseNameA, caseExtA, smooth, window)
dataSetE, metaDataE = dataExtraction(workingDirE, caseNameA, caseExtA, smooth, window)
metaData = dataSetPlot(dataSetA, metaData, dA, smooth=window)
# metaData = dataSetPlot(dataSetB, metaData, dB, smooth=window)
# metaData = dataSetPlot(dataSetC, metaData, dC, smooth=window)
# metaData = dataSetPlot(dataSetD, metaData, dD, smooth=window)
metaData = dataSetPlot(dataSetE, metaData, dE, smooth=window)
metaPlot(metaData, prop='ReP')


ax1.set_title("PDFs")
ax2.set_title("PDFs")
ax3.set_title("Differentiated PDF")
ax1.set_xlabel("Value")
ax1.set_ylabel("Normalized freq.")
ax2.set_xlabel("Value")
ax2.set_ylabel("Normalized freq.")
ax1.legend(loc=0)
ax1.set_yscale('log')
ax2.legend(loc=0)
ax3.legend(loc=0)

sns.despine(f1)
sns.despine(f2)
sns.despine(f3)
sns.despine(f4)
sns.despine(f5)
sns.despine(f6)

# ax3.set_yscale('log')
# plt.xscale('log')
plt.ion()
plt.show()
