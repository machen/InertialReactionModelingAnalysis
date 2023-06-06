import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
import seaborn as sns
from itertools import cycle

"""
TODO: Data that is collated should be output as a csv (or at least have a flag to output as a CSV)
TODO: Implement dataset method which would allow for sweeping an arbitrary number of folders

TODO: I want to look at variation in images which means bringing in data images that have the same conditions.
np.stack allows me to do this by creating a 3D array. I can then get the "mean image" and "std image" from that array and output it.

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




def flowRateConversion(q, width, height, charLen, nu=0.43):
    # Base nu is in mm2/s, so you should report this in mm and seconds as units
    reyn = q/width/height*charLen/nu
    return reyn


def extractExpParams(fileName):
    # Produces a dictionary of experimental parameters
    qPat = re.compile('(\d*\.?\d*)q')
    repPat = re.compile('(\d{0,3}).nd2')
    qVal = re.search(qPat, fileName).group(1)
    repVal = re.search(repPat, fileName).group(1)
    if not repVal:
        repVal = 0
    res = {'q': float(qVal), 'replicate': int(repVal)}
    return res


def loadMetaData(metaFile, caseExt):
    params = pd.read_csv(metaFile, index_col='fileName')
    params = params.loc[params.caseExt == caseExt]
    return params


def dataExtraction(workingDir, caseName, caseExt, smooth=False, window=5):
    filePat = re.compile(caseName+r".*"+caseExt)
    fileList = os.listdir(workingDir)
    dataSets = {}
    for fileName in fileList:
        if re.match(filePat, fileName):
            print(fileName)
            # params = extractExpParams(fileName)
            data = pd.read_csv(workingDir+fileName, header=0)
            if smooth:
                # Smooth data with rolling average of size window
                data = dataSmoothing(data, window)
                # Drop NaNs
                data = data.dropna()
            dataSets[os.path.splitext(fileName)[0]] = data
    metaData = loadMetaData(workingDir+'_meta.csv', caseExt)
    return dataSets, metaData


def dataSetPlot(dataSets, metaData, d, linestyle='-', smooth=0):
    for key in dataSets:
        data = dataSets[key]
        dataMean, dataVar = pdfStats(data)
        params = extractExpParams(key)
        # params['d'] = d
        # params['fileName'] = key
        # params['PDFmean'] = dataMean
        # params['PDFstd'] = np.sqrt(dataVar)
        metaData.loc[key, 'd'] = d
        metaData.loc[key, 'PDFmean'] = dataMean
        metaData.loc[key, 'PDFstd'] = np.sqrt(dataVar)
        #metaData = metaData.append(params, ignore_index=True)
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


def metaPlot(metaData, prop='q', marker='o', propLim=None):
    # Change me to work with the image histograms and metadata
    if prop == 'Re':
        metaData.loc[:, prop] = flowRateConversion(metaData.q/60.0, 250E-3, 100E-3, 250E-3)
    elif prop == 'ReP':
        metaData.loc[:, prop] = flowRateConversion(metaData.q/60.0, 250E-3, 100E-3, 100E-3)
    metaData.sort_values(by=prop)
    traceList = metaData.d.unique()
    for val in traceList:
        subData = metaData.loc[metaData.d == val, :]
        if propLim:
            subData = subData.loc[subData.loc[:, prop] < propLim, :]
        meanData = subData.groupby(prop).mean()
        stdData = subData.groupby(prop).std()
        ax4.plot(subData.loc[:, prop], subData.sumInt, ls='none', marker=marker,
                 label=val)
        ax5.plot(subData.loc[:, prop], subData.loc[:, 'reactorRatio'],
                 ls='none', marker=marker, label=val)
        maxLoc = meanData.loc[meanData.meanInt ==
                              meanData.meanInt.max(), 'q'].values
        ax6.errorbar(meanData.index, meanData.loc[:, 'meanInt'],
                     yerr=stdData.loc[:,'stdInt'], ls='none', marker=marker,
                     capsize=2, label="{} max val q: {}".format(val,
                                                                float(maxLoc)))
        maxVal = meanData.PDFmean.max()
        ax7.plot(meanData.index, meanData.loc[:, 'reactorRatio'],
                 ls='none', marker=marker, label=val)
        maxVal = meanData.sumInt.max()
        ax8.errorbar(meanData.index, meanData.sumInt/maxVal,
                     yerr=stdData.sumInt/maxVal, ls='-',
                     marker=marker, capsize=2, label=val)
        ax9.errorbar(meanData.index, meanData.sumInt, yerr=stdData.sumInt, ls='none',
                     marker=marker, capsize=2, label=val)
    ax4.set_xlabel(prop)
    ax4.set_ylabel('Sum intensity individual image')
    ax4.legend(loc=0)
    ax5.set_xlabel(prop)
    ax5.set_ylabel('Reactor ratio')
    ax5.legend(loc=0)
    ax6.set_xlabel(prop)
    ax6.set_ylabel('Mean Int')
    ax6.set_title('Averaged Results')
    ax6.legend(loc=0)
    ax7.set_xlabel(prop)
    ax7.set_ylabel('Reactor ratio, averaged over images')
    ax7.legend(loc=0)
    ax8.legend(loc=0)
    ax8.set_xlabel(prop)
    ax8.set_ylabel('Sum intensity averaged over images')
    ax8.set_title('Max normalized - check meaning of max val')
    ax8.set_ylim([0.3, 1.05])
    # ax8.set_xlim([-1, 105])
    ax9.legend(loc=0)
    ax9.set_xlabel(prop)
    ax9.set_ylabel('Sum Intensity averaged over images')
    ax9.set_title('Sum intensities')
    # ax9.set_xlim([-1, 105])
    meanDict = {'q':meanData.q.values,'meanInt':meanData.meanInt.values,'stdMeanInt':stdData.meanInt.values,
           'normMeanInt':meanData.meanInt.values/maxVal, 'stdNormMeanInt':stdData.meanInt.values/maxVal}
    out = pd.DataFrame(meanDict)
    out.loc[:,'ID'] = val
    return out
sns.set_context('poster', font_scale=1.25)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Cambria'

smooth = False
window = 10
prop = 'ReP'
propLim = None
# Might be nice to do some averaging of lines that have the same experiemntal condition

mainDir = "..\\..\\Experiments\\"

# workingDirA = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Image 1 50 Bins\\"
# workingDirB = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Image 2 50 Bins\\"
# workingDirC = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Image 3 50 Bins\\"
# workingDirD = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Image 4 50 Bins\\"
# workingDirE = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Whole Channel Aligned 50 Bins\\"
# workingDirF = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 7 50 Bins\\"
# workingDirG = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 8 50 Bins\\"
# workingDirH = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 9 50 Bins\\"
# workingDirI = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 10 50 Bins\\"
# workingDirJ = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 11 50 Bins\\"
# workingDirK = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 12 50 Bins\\"
# workingDirL = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 13 50 Bins\\"
# workingDirM = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Pore 14 50 Bins\\"
# workingDirN = "2023-3-21-Chemilum\\MPD1_D4_Batch2\\Whole Channel Aligned 50 Bins\\"

# workingDirA = "2023-2-25 Chemilum\\MPD1_D4\\Image 1 Aligned Images 50 Bins\\"
# workingDirB = "2023-2-25 Chemilum\\MPD1_D4\\Image 2 Aligned Images 50 Bins\\"
# workingDirC = "2023-2-25 Chemilum\\MPD1_D4\\Image 3 Aligned Images 50 Bins\\"
# workingDirD = "2023-2-25 Chemilum\\MPD1_D4\\Image 4 Aligned Images 50 Bins\\"
# workingDirE = "2023-2-25 Chemilum\\MPD1_D4\\Stitched Aligned Images 50 Bins\\"

workingDirA = "2023-5-26-Tracer_Chemilum\\Tracer\\Top half Conc 50 bins\\"
workingDirB = "2023-5-26-Tracer_Chemilum\\Tracer\\Bottom half Conc 50 bins\\"
workingDirC = "2023-5-26-Tracer_Chemilum\\Tracer\\Whole region Conc 50 bins\\"
# workingDirG = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 6 50 bins\\"
# workingDirH = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 7 50 bins\\"
# workingDirI = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 8 50 bins\\"
# workingDirJ = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 9 50 bins\\"
# workingDirK = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 10 50 bins\\"
# workingDirL = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 11 50 bins\\"
# workingDirM = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 12 50 bins\\"
# workingDirN = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 13 50 bins\\"
# workingDirO = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 14 50 bins\\"
# workingDirP = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 15 50 bins\\"
# workingDirQ = "2022-3-22-MPD2\\MPD2_P1_A3\\S1 Raw Masked Pore Throat 16 50 bins\\"

dA = "Top half"
dB = "Bottom half"
dC = "Whole channel"
dD = "Batch1 Mid Pore B"
dE = "Batch1 Inlet Pore C"
dF = "Batch 1 Inlet Pores"
# dG = "Pore 6"
# dH = "Pore 7"
# dI = "Pore 8"
# dJ = "Pore 9"
# dK = "Pore 10"
# dL = "Pore 11"
# dM = "Pore 12"
# dN = "Pore 13"
# dO = "Pore 14"
# dP = "Pore 15"
# dQ = "Pore 16"


# workingDirA = "2022-5-19-Chemilum\\2PD4_P7_A2\\Pillar Gap Exclusive 50 bins\\"
# workingDirB = "2022-1-15-Chemilum 100 um\\100 um Gap\\Pillar Gap Exclusive Aligned 50 bins\\"
# workingDirC = "2021-10-20-Chemilum-100um\\A3-100um\\Pillar Gap Exclusive Aligned 50 bins\\"
# workingDirD = "2022-2-1-Chemilum\\2PD1_P7_A2\\Pillar Gap Exclusive Aligned 50 bins\\"
# workingDirE = "2021-10-05-Chemilum-100um\\100 um Pillar Gap\\Pillar Gap Exclusive Aligned 50 bins\\"
# workingDirF = "2021-11-18-Chemilum-25um\\2PD3_A2\\Pillar Gap Exclusive Aligned 50 bins\\"
# workingDirG = "2022-2-9-Chemilum\\25umGap\\Pillar Gap Exclusive 50 bins\\"
# dA = "1/15"
# dB = "10/5"
# dC = "10/20"
# dD = "2/1"

# workingDirA = "2022-2-9-Chemilum\\Multipillar\\Raw Aligned Image Pillar Gap Image 1 50 bins\\"
# workingDirB = "2022-2-9-Chemilum\\Multipillar\\Raw Aligned Image Pillar Gap Image 3 50 bins\\"
# workingDirC = "2022-2-9-Chemilum\\Multipillar\\Raw Aligned Image Pillar Gap Image 5 50 bins\\"

# workingDirA = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Center Pore Throats 50 bins\\"
# workingDirB = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore 3 50 bins\\"
# workingDirC = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore 9 50 bins\\"
# workingDirD = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore 10 50 bins\\"
# workingDirE = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore 5 50 bins\\"
# workingDirF = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore 9 50 bins\\"
# workingDirG = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore 15 50 bins\\"
# # workingDirC = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Whole Channel 50 bins\\"
# workingDirA = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore Throats 50 bins\\"
# workingDirD = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Pore 3 50 bins\\"

# workingDirA = "2022-3-22-MPD2\\MPD2_P1_A3\\Raw Masked Image 50 bins\\"


# dF = "Image 5"
# dG = "2/9-25 um"
# dF = "Pore 9"
# dG = "Pore 15"
# dE = "10/20 A2"

# workingDirA = "2022-1-15-Chemilum 100 um\\100 um Gap\\Raw Aligned Image Pillar Gap 50 bins\\"
# workingDirB = "2022-2-9-Chemilum\\25umGap\\Raw Aligned Image Pillar Gap 50 bins\\"
# workingDirC = "2021-10-20-Chemilum-100um\\A3-100um\\Raw Aligned Image Pillar Gap 50 bins\\"
# workingDirD = "2022-2-1-Chemilum\\2PD1_P7_A2\\Raw Aligned Image Pillar Gap 50 bins\\"
# workingDirA = "2021-11-18-Chemilum-25um\\2PD3_A2\\Raw Aligned Image Pillar Gap 50 bins\\"
# workingDirB = "2021-10-05-Chemilum-100um\\100 um Pillar Gap\\Raw Aligned Image Pillar Gap 50 bins\\"
# workingDirC = "2021-05-19-Chemilum-25um\\ExptImages\\Raw Aligned Image Pillar Gap 50 bins\\"
# workingDirH = "2022-1-6-Chemilum\\100 um Gap\\Raw Aligned Image Pillar Gap 50 bins\\"

os.chdir(mainDir)
caseNameA = ''
caseExtA = r".fluor_hist" # TODO: YOU NEED TO DROP METADATA BASED ON WHICH CHANNEL YOU ARE SELECTING

# You must set these to the correct pillar gaps of the experiment

# dF = "3/22/2022 Image 5"
# dG = "3/22/2022 Whole Channel"
# dH = "1/6 100 um"

# dD = "2022-2-1 100 um"
# dE = "2022-10-20 100 um"

f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f3, ax3 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f4, ax4 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f5, ax5 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f6, ax6 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f7, ax7 = plt.subplots(1, 1, sharex='col', figsize=(12, 10)) # Plot for normalizing to max average value
f8, ax8 = plt.subplots(1, 1, sharex='col', figsize=(12, 10)) # Plot using meanIntensity rather than mean of the PDF
f9, ax9 = plt.subplots(1, 1, sharex='col', figsize=(12, 10)) # Plot using meanIntensity rather than mean of the PDF
f10, ax10 = plt.subplots(1,1, sharex='col', figsize=(12, 10))

metaData = pd.DataFrame([], columns=['q', 'replicate', 'PDFmean', 'PDFstd'])
dataSetA, metaDataA = dataExtraction(workingDirA, caseNameA, caseExtA, smooth, window)
dataSetB, metaDataB = dataExtraction(workingDirB, caseNameA, caseExtA, smooth, window)
dataSetC, metaDataC = dataExtraction(workingDirC, caseNameA, caseExtA, smooth, window)
# dataSetD, metaDataD = dataExtraction(workingDirD, caseNameA, caseExtA, smooth, window)
# dataSetE, metaDataE = dataExtraction(workingDirE, caseNameA, caseExtA, smooth, window)
# dataSetF, metaDataF = dataExtraction(workingDirF, caseNameA, caseExtA, smooth, window)
# dataSetG, metaDataG = dataExtraction(workingDirG, caseNameA, caseExtA, smooth, window)
# dataSetH, metaDataH = dataExtraction(workingDirH, caseNameA, caseExtA, smooth, window)
# dataSetI, metaDataI = dataExtraction(workingDirI, caseNameA, caseExtA, smooth, window)
# dataSetJ, metaDataJ = dataExtraction(workingDirJ, caseNameA, caseExtA, smooth, window)
# dataSetK, metaDataK = dataExtraction(workingDirK, caseNameA, caseExtA, smooth, window)
# dataSetL, metaDataL = dataExtraction(workingDirL, caseNameA, caseExtA, smooth, window)
# dataSetM, metaDataM = dataExtraction(workingDirM, caseNameA, caseExtA, smooth, window)
# dataSetN, metaDataN = dataExtraction(workingDirN, caseNameA, caseExtA, smooth, window)
# dataSetO, metaDataO = dataExtraction(workingDirO, caseNameA, caseExtA, smooth, window)
# dataSetP, metaDataP = dataExtraction(workingDirP, caseNameA, caseExtA, smooth, window)
# dataSetQ, metaDataQ = dataExtraction(workingDirQ, caseNameA, caseExtA, smooth, window)


metaDataA = dataSetPlot(dataSetA, metaDataA, dA, smooth=window)
metaDataB = dataSetPlot(dataSetB, metaDataB, dB, smooth=window)
metaDataC = dataSetPlot(dataSetC, metaDataC, dC, smooth=window)
# metaDataD = dataSetPlot(dataSetD, metaDataD, dD, smooth=window)
# metaDataE = dataSetPlot(dataSetE, metaDataE, dE, smooth=window)
# metaDataF = dataSetPlot(dataSetF, metaDataF, dF, smooth=window)
# metaDataG = dataSetPlot(dataSetG, metaDataG, dG, smooth=window)
# metaDataH = dataSetPlot(dataSetH, metaDataH, dH, smooth=window)
# metaDataI = dataSetPlot(dataSetI, metaDataI, dI, smooth=window)
# metaDataJ = dataSetPlot(dataSetJ, metaDataJ, dJ, smooth=window)
# metaDataK = dataSetPlot(dataSetK, metaDataK, dK, smooth=window)
# metaDataL = dataSetPlot(dataSetL, metaDataL, dL, smooth=window)
# metaDataM = dataSetPlot(dataSetM, metaDataM, dM, smooth=window)
# metaDataN = dataSetPlot(dataSetN, metaDataN, dN, smooth=window)
# metaDataO = dataSetPlot(dataSetO, metaDataO, dO, smooth=window)
# metaDataP = dataSetPlot(dataSetP, metaDataP, dP, smooth=window)
# metaDataQ = dataSetPlot(dataSetQ, metaDataQ, dQ, smooth=window)

markerCycle = cycle(['o', 'd', 's', '^', 'D', 'h', 'X'])
out = metaPlot(metaDataA, prop=prop, marker=next(markerCycle), propLim=propLim)
out = pd.concat([out, metaPlot(metaDataB, prop=prop, marker=next(markerCycle), propLim=propLim)])
out = pd.concat([out, metaPlot(metaDataC, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataD, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataE, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataF, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataG, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataH, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataI, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataJ, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataK, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataL, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataM, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataN, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataO, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataP, prop=prop, marker=next(markerCycle), propLim=propLim)])
# out = pd.concat([out, metaPlot(metaDataQ, prop=prop, marker=next(markerCycle), propLim=propLim)])




# metaPlot(metaDataG, prop=prop, marker=next(markerCycle))
# metaPlot(metaDataH, prop=prop, marker=next(markerCycle))

out = out.dropna()
out.loc[:,'ReP'] = flowRateConversion(out.q/60.0, 250E-3, 100E-3, 100E-3)
traceAvg = out.dropna().groupby('q').mean()
traceStd = out.dropna().groupby('q').std()
ax10.errorbar(traceAvg.ReP, traceAvg.normMeanInt, traceStd.normMeanInt, color='k', marker='o')
ax10.set_title('Averaged Traces')
ax10.set_xlabel('Reynolds Number')
ax10.set_ylabel('Normalized Intensity (.)')
ax10.set_ylim([0.3, 1.05])
ax10.set_xlim([-1, 105])

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
sns.despine(f7)
sns.despine(f8)
sns.despine(f9)

# ax3.set_yscale('log')
# plt.xscale('log')
plt.ion()
plt.show()
