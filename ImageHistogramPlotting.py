import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
# import seaborn as sns


def extractParams(fileName):
    # Produces a dictionary of experimental parameters
    qPat = re.compile('(\d*\.?\d*)q')
    repPat = re.compile('(\d{0,3})_dark_hist.csv')
    qVal = re.search(qPat, fileName)
    repVal = re.search(repPat, fileName)
    res = {'q': float(qVal.group(1)), 'replicate': repVal.group(1)}
    return res


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
    return dataSets


def dataSetPlot(dataSets, metaData, linestyle='-', smooth=0, fit=True):
    for key in dataSets:
        data = dataSets[key]
        params = extractParams(key)
        dataMean, dataVar = pdfStats(data)
        params['PDFmean'] = dataMean
        params['PDFstd'] = np.sqrt(dataVar)
        metaData = metaData.append(params, ignore_index=True)
        a1 = ax1.plot(data.valMean, data.normFreq,
                      label=key+'smooth {}'.format(smooth), ls=linestyle)
        c = a1[0].get_color()
        ax1.plot([dataMean, dataMean],
                 [np.min(data.normFreq), np.max(data.normFreq)],
                 ls='-', color='k')
        # ax1.plot([dataMean-np.sqrt(dataVar), dataMean+np.sqrt(dataVar)],
        #          [dataMid, dataMid], ls='--', color=c)

        a2 = ax2.plot(data.valMean,
                      data.normFreq, label=key+'smooth {}'.format(smooth),
                      ls=linestyle)
        c = a2[0].get_color()
        ax2.plot([dataMean, dataMean],
                 [np.min(data.normFreq), np.max(data.normFreq)],
                 ls='-', color='k')
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
    metaData.sort_values(by=prop)
    ax4.errorbar(metaData.loc[:, prop], metaData.loc[:, 'PDFmean'],
                 yerr=metaData.loc[:, 'PDFstd'], ls='none', marker='o',
                 capsize=2)
    return


plt.rcParams['svg.fonttype'] = 'none'
smooth = False
window = 10


workingDirA = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\2PillarD-1_P2_A2\\Result 100 bins\\"
os.chdir(workingDirA)
caseNameA = '\d*q_Pillars\d*'

caseExtA = ".*dark_hist\.csv"

f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f3, ax3 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f4, ax4 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))

metaData = pd.DataFrame([], columns=['q', 'replicate', 'PDFmean', 'PDFstd'])
dataSetA = dataExtraction(workingDirA, caseNameA, caseExtA, smooth, window)
metaData = dataSetPlot(dataSetA, metaData, smooth=window)
metaPlot(metaData)
# dataSetB = dataExtraction(workingDirB, caseNameB, caseExtB, smooth, window)
# metaData = dataSetPlot(dataSetB, metaData, smooth=window)

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
ax4.set_xlabel('q')
ax4.set_ylabel('Mean of PDF')
ax3.set_yscale('log')
# plt.xscale('log')
ax4.legend(loc=0)
plt.ion()
plt.show()