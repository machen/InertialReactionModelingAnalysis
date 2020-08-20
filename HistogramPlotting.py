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


def dataPlot(workingDir, caseName, caseExt):
    os.chdir(workingDir)
    filePat = re.compile(caseName+'.*?'+caseExt)
    fileList = os.listdir('.')
    for fileName in fileList:
        if re.match(filePat, fileName):
            print(fileName)
            # params = extractParams(fileName)
            data = pd.read_csv(fileName, header=0,
                               names=['binID', 'normFreq', 'valMean'])
            ax1.plot(data.valMean, data.normFreq, label=fileName)
            ax2.plot(data.valMean/data.valMean.mean(),
                     data.normFreq, label=fileName)
    return


# workingDirA = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\ChemData\\PillarGap_Norm"
workingDirA = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData_FlowOnly\\Pillar gap - 500 log bins\\"
# workingDir = "."
caseNameA = "TwoInletsTwoColumns_v5.2_ExF_FlowOnly"
caseExtA = "Re250.flowdata_histogram\.csv$"

workingDirB = "..\\..\\..\\..\\..\\Multipillar\\Normal\\FlowData_Normal\\200 log bins - 250 to -2500"
caseNameB = "Multipillar_v5.2_Normal_r100"
caseExtB = "Re100_histogram\.csv$"

# Plot for everything
f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))

dataPlot(workingDirA, caseNameA, caseExtA)
# dataPlot(workingDirB, caseNameB, caseExtB)


ax1.set_title("PDFs")
ax2.set_title("Normalized velocity PDFs")
ax1.set_xlabel("Value")
ax1.set_ylabel("Normalized freq.")
ax2.set_xlabel("Val/Average Val")
ax2.set_ylabel("Normalized freq.")
ax1.legend(loc=0)
ax2.legend(loc=0)
plt.ion()
plt.show()
