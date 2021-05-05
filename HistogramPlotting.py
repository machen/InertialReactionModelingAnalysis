import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
import seaborn as sns

"""Purpose of script is to open and plot generated histogram files, skipping
the analysis, but instead just putting out the images.


CURRENTLY SWEEPS A WORKING FOLDER, TURN INTO A FUNCTION

WOULD BE NICE TO MAKE A COMPARITOR FUNCTION THAT CAN TAKE TWO FILE NAMES (OR TWO COMPARISON DIRECTORIES)


HEY SOME STUFF YOU NEED TO DO

-Calculate distribution median (point at which each area is divided half and half)
-create a plot that compares various aggergate statistics to a given parameter


 """
plt.rcParams['svg.fonttype'] = 'none'
sns.set_context('talk')


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
        r1Pat = re.compile('r1_(\d+?)_')
        r2Pat = re.compile('r2_(\d+?)_')
        rePat = re.compile('Re(.*?)\.(chem|flow)data')
        dPat = re.compile('d(\d+?)_')
        cPat = re.compile('c(\d+?)')
        kPat = re.compile('k(\d+?)_')
        r1Val = float(re.search(r1Pat, fileName).group(1))
        r2Val = float(re.search(r2Pat, fileName).group(1))
        dVal = float(re.search(dPat, fileName).group(1))
        reVal = float(re.search(rePat, fileName).group(1))
        cVal = float(re.search(cPat, fileName).group(1))
        kVal = float(re.search(kPat, fileName).group(1))
        vel = reVal*1.6E-6/500E-6  # Viscosity, and 500 um channel width
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
    rePillar = res['Re']*res['r1']*2/500  # Assume 500 um channel width
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


def dataSetPlot(dataSets, metaData, linestyle='-', smooth=0, fit=True):
    # In this case, key refers to the filename for the given data set
    workingDirParams = importParams()
    # print(workingDirParams)
    for key in dataSets:
        print(key)
        data = dataSets[key]
        params = extractParams(key, workingDirParams)
        dataMean, dataVar = pdfStats(data)
        params['PDFmean'] = dataMean
        params['PDFstd'] = np.sqrt(dataVar)
        metaData = metaData.append(params, ignore_index=True)
        dataMid = np.mean([np.min(data.normFreq), np.max(data.normFreq)])
        pmf, xPmf = genPMF(data)
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
        ax5.plot(xPmf, pmf, label=key, ls=linestyle)
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


def metaPlot(metaData, prop='Re', flowCond='NS'):
    subData = metaData.loc[metaData.loc[:, 'Flow'] == flowCond, :]
    subData.sort_values(by=prop)
    ax4.errorbar(subData.loc[:, prop], subData.loc[:, 'PDFmean'],
                 yerr=subData.loc[:, 'PDFstd'], ls='none', marker='o',
                 capsize=2, label=flowCond+" "+prop)
    ax4.set_xlabel(prop)
    ax6.errorbar(subData.loc[:, prop], subData.loc[:, 'RePil'],
                 yerr=subData.loc[:, 'PDFstd'], ls='none', marker='o',
                 capsize=2, label=flowCond+" "+prop)
    ax6.set_xlabel(prop)
    ax6.set_ylabel('Re')
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


quickVol = {0.1: 1.48E-13, 1: 6.02E-14, 10: 5.85E-14, 50: 4.58E-13}

window = 5
smooth = False
fitRange = np.array([85, 90])
prop = 'RePil'  # Options: Re, d, RePil, DaAdv, DaDiff, Pe, reactorConserv
# fitRange = np.array([65, 85])
workingDirA = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\FlowData\\Pillar Gap Exact-angle-180 linear bins\\"
#workingDirA = "..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\Pillar Gap-angle-180 linear bins"
#workingDirA = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\ChemData\\Pillar Gap-dCdtNorm-100 linear bins\\"
# workingDir = "."
caseNameA = "TwoPillar_v6"
caseExtA = "d100_Re.*\.flowdata_histogram\.csv"
# workingDirB = "..\\..\\..\\..\\..\\Multipillar\\Normal\\FlowData_Normal\\200 log bins - 250 to -2500"
#workingDirB = "..\\..\\..\\..\\..\\..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\Pillar gap-angle-180 linear bins"
workingDirB = "..\\..\\..\\..\\..\\..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\FlowData\\Pillar Gap Exact-angle-180 linear bins\\"
caseNameB = "TwoPillar_v6"
caseExtB = "THING\.flowdata_histogram\.csv"

# Plot for everything
f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f3, ax3 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f4, ax4 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f5, ax5 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f6, ax6 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))

metaData = pd.DataFrame([], columns=['r1', 'r2', 'd', 'Re', 'Flow', 'PDFmean', 'PDFstd'])
dataSetA = dataExtraction(workingDirA, caseNameA, caseExtA, smooth, window)
metaData = dataSetPlot(dataSetA, metaData, smooth=window, linestyle='-')
# for key in dataSetA:
#     fit, xVal, yRes = semilogYFitting(dataSetA[key], 'valMean',
#                                       'normFreq', np.array([87, 91]))
#     ax1.plot(xVal, yRes, ls='-', color='k', label='A'+str(fit))
#     ax2.plot(xVal, yRes, ls='-', color='k', label='A'+str(fit))
#     fit, xVal, yRes = semilogYFitting(dataSetA[key], 'valMean',
#                                       'normFreq', np.array([75, 87]))
#     ax1.plot(xVal, yRes, ls='-', color='k', label='A'+str(fit))
#     ax2.plot(xVal, yRes, ls='-', color='k', label='A'+str(fit))
# for i in [2,4,5,10]:
#     dataSetSmooth = dataExtraction('.', caseNameA, caseExtA, smooth=True, window=i)
#     dataSetPlot(dataSetSmooth, smooth=i)
dataSetB = dataExtraction(workingDirB, caseNameB, caseExtB, smooth, window)
metaData = dataSetPlot(dataSetB, metaData, smooth=window,linestyle='--')
# for key in dataSetB:
#     fit, xVal, yRes = semilogYFitting(dataSetB[key], 'valMean',
#                                       'normFreq', np.array([75, 91]))
#     ax1.plot(xVal, yRes, ls='--', color='k', label='B'+str(fit))
#     ax2.plot(xVal, yRes, ls='--', color='k', label='B'+str(fit))
# for i in range(10):
#     xGauss, yGauss = genGaussian(90, i)
#     ax1.plot(xGauss, yGauss, ls='--', color='k', label='{} gaussian'.format(i))
#     ax2.plot(xGauss, yGauss, ls='--', color='k', label='{} gaussian'.format(i))

metaPlot(metaData, prop=prop, flowCond='NS')
metaPlot(metaData, prop=prop, flowCond='Stokes')

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

ax4.set_ylabel('Mean of PDF')
ax3.set_yscale('log')
# plt.xscale('log')
ax4.legend(loc=0)
ax5.set_title("PMFs")
ax5.set_xlabel('Value')
ax5.set_ylabel('PMF')
ax5.legend(loc=0)
ax5.set_yscale('log')

sns.despine(f1)
sns.despine(f2)
sns.despine(f3)
sns.despine(f4)
sns.despine(f5)

plt.ion()
plt.show()
