import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
# import seaborn as sns

"""Purpose of script is to open and plot generated histogram files, skipping
the analysis, but instead just putting out the images.


CURRENTLY SWEEPS A WORKING FOLDER, TURN INTO A FUNCTION

WOULD BE NICE TO MAKE A COMPARITOR FUNCTION THAT CAN TAKE TWO FILE NAMES (OR TWO COMPARISON DIRECTORIES)

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
            data = pd.read_csv(fileName, header=0,
                               names=['binID', 'normFreq', 'valMean'])
            if smooth:
                # Smooth data with rolling average of size window
                data = dataSmoothing(data, window)
                # Drop NaNs
                data = data.dropna()
            dataSets[fileName] = data
    return dataSets


def dataSetPlot(dataSets, linestyle='-'):
    for key in dataSets:
        data = dataSets[key]
        ax1.plot(data.valMean, data.normFreq, label=key, ls=linestyle)
        ax2.plot(data.valMean,
                 data.normFreq, label=key, ls=linestyle)
    return


def dataSmoothing(data, window=5):
    dataSmooth = data.rolling(window, center=True).mean()
    return dataSmooth

window = 5
smooth = True

workingDirA = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData\\Pillar region-angle-180 linear bins\\"
#workingDirA = "..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\Pillar Region - Angle - 180 linear bins"
# workingDir = "."
caseNameA = "TwoInletsTwoColumns_v5.2_"
caseExtA = "Re10.flowdata_histogram\.csv"

# workingDirB = "..\\..\\..\\..\\..\\Multipillar\\Normal\\FlowData_Normal\\200 log bins - 250 to -2500"
workingDirB = "..\\..\\..\\..\\..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData_FlowOnly\\Pillar region - angle - 180 linear bins"
caseNameB = "TwoInletsTwoColumns_v5.2_ExF_FlowOnly_GapVar_r1_100_r2_100_d50"
caseExtB = ".flowdata_histogram\.csv"

# Plot for everything
f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))

dataSetA = dataExtraction(workingDirA, caseNameA, caseExtA, smooth, window)
dataSetPlot(dataSetA)
# dataSetB = dataExtraction(workingDirB, caseNameB, caseExtB)
# dataPlot(dataSetB)
# ax2.plot([0, 0], [0.01, 1E6], ls='--', color='k')

ax1.set_title("PDFs")
ax2.set_title("PDFs")
ax1.set_xlabel("Value")
ax1.set_ylabel("Normalized freq.")
ax2.set_xlabel("Value")
ax2.set_ylabel("Normalized freq.")
ax1.legend(loc=0)
ax1.set_yscale('log')
ax2.legend(loc=0)
plt.yscale('log')
plt.xscale('log')
plt.ion()
plt.show()
