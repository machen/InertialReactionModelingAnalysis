import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import scipy.stats as stats
# import seaborn as sns

"""Purpose of script is to open and plot generated histogram files, skipping
the analysis, but instead just putting out the images.


CURRENTLY SWEEPS A WORKING FOLDER, TURN INTO A FUNCTION

WOULD BE NICE TO MAKE A COMPARITOR FUNCTION THAT CAN TAKE TWO FILE NAMES (OR TWO COMPARISON DIRECTORIES)


HEY SOME STUFF YOU NEED TO DO
-Estimate a derivative of a given data set (np.diff(freq)/np.diff(valMean))
-Calculate the mean of the distribution (integral of valMean*freq dval)
-Calculate distribution median (point at which each area is divided half and half)

 """
plt.rcParams['svg.fonttype'] = 'none'


def extractParams(fileName):
    # Produces a dictionary of experimental parameters
    r1Pat = re.compile('r1_(\d+?)_')
    r2Pat = re.compile('r2_(\d+?)_')
    rePat = re.compile('Re(.*?)_')
    dPat = re.compile('d(\d+?)_')
    r1Val = re.search(r1Pat, fileName).group(1)
    r2Val = re.search(r2Pat, fileName).group(1)
    dVal = re.search(dPat, fileName).group(1)
    reVal = re.search(rePat, fileName).group(1)
    res = {'r1': float(r1Val), 'r2': float(r2Val), 'd': float(dVal),
           'Re': float(reVal)}
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


def dataSetPlot(dataSets, linestyle='-', smooth=0, fit=True):
    for key in dataSets:
        data = dataSets[key]
        ax1.plot(data.valMean, data.normFreq,
                 label=key+'smooth {}'.format(smooth), ls=linestyle)
        ax2.plot(data.valMean,
                 data.normFreq, label=key+'smooth {}'.format(smooth),
                 ls=linestyle)
# HEY UNIFY THE NAMING IN THE MAIN SCRIPT SINCE IT'S NOT ALWAYS VELOCITIES
    return


def dataSmoothing(data, window=5):
    dataSmooth = data.rolling(window, center=True).mean()
    return dataSmooth


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

window = 5
smooth = True
fitRange = np.array([85, 90])
#fitRange = np.array([65, 85])

workingDirA = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData\\Pillar Gap-angle-180 linear bins\\"
#workingDirA = "..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\Pillar Gap -angle- 180 linear bins"
#workingDirA = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\ChemData\\Pillar gap-cProduct-100 linear bins\\"
# workingDir = "."
caseNameA = "TwoInletsTwoColumns_v5.2_ExF_r100_d100"
caseExtA = ".chemdata_histogram\.csv"
# workingDirB = "..\\..\\..\\..\\..\\Multipillar\\Normal\\FlowData_Normal\\200 log bins - 250 to -2500"
workingDirB = "..\\..\\..\\..\\..\\..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\Pillar gap-angle-180 linear bins"
#workingDirB = "..\\..\\..\\..\\..\\..\\Comsol5.5\\TwoPillars\\ExF\\ChemData\\Pillar gap-cProduct-100 linear bins\\"
caseNameB = "TwoInletsTwoColumns_v5.2_ExF_r100_d100"
caseExtB = ".chemdata_histogram\.csv"

# Plot for everything
f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))


dataSetA = dataExtraction(workingDirA, caseNameA, caseExtA, smooth, window)
dataSetPlot(dataSetA, smooth=window)
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
dataSetPlot(dataSetB, smooth=window)
# for key in dataSetB:
#     fit, xVal, yRes = semilogYFitting(dataSetB[key], 'valMean',
#                                       'normFreq', np.array([75, 91]))
#     ax1.plot(xVal, yRes, ls='--', color='k', label='B'+str(fit))
#     ax2.plot(xVal, yRes, ls='--', color='k', label='B'+str(fit))
# for i in range(10):
#     xGauss, yGauss = genGaussian(90, i)
#     ax1.plot(xGauss, yGauss, ls='--', color='k', label='{} gaussian'.format(i))
#     ax2.plot(xGauss, yGauss, ls='--', color='k', label='{} gaussian'.format(i))

ax1.set_title("PDFs")
ax2.set_title("PDFs")
ax1.set_xlabel("Value")
ax1.set_ylabel("Normalized freq.")
ax2.set_xlabel("Value")
ax2.set_ylabel("Normalized freq.")
ax1.legend(loc=0)
ax1.set_yscale('log')
ax2.legend(loc=0)
# plt.yscale('log')
# plt.xscale('log')
plt.ion()
plt.show()
