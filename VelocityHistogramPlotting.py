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
    rePat = re.compile('Re(.*?)_histogram.csv')
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
                               names=['binID', 'normFreq', 'velVal'])
            ax1.plot(data.velVal, data.normFreq, label=fileName)
            ax2.plot(data.velVal/data.velVal.mean(),
                     data.normFreq, label=fileName)
    return


workingDirA = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData_FlowOnly\\Pillar region - 500 log bins"
# workingDir = "."
caseNameA = "TwoInletsTwoColumns_v5.1_ExF_FlowOnly_r1_100"
caseExtA = "100_histogram\.csv$"

workingDirB = "..\\..\\..\\..\\..\\Multipillar\\Normal\\FlowData_Normal\\200 log bins - 250 to -2500"
caseNameB = "Multipillar_v5.2_Normal_r100"
caseExtB = "Re100_histogram\.csv$"


# Plot for everything
f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
f2, ax2 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))

dataPlot(workingDirA, caseNameA, caseExtA)
dataPlot(workingDirB, caseNameB, caseExtB)


ax1.set_title("Velocity PDFs")
ax2.set_title("Normalized velocity PDFs")
ax1.set_xlabel("Velocity (m/s)")
ax1.set_ylabel("Normalized freq.")
ax2.set_xlabel("Velocity/Average Velocity")
ax2.set_ylabel("Normalized freq.")
ax1.legend(loc=0)
ax2.legend(loc=0)
plt.ion()
plt.show()
