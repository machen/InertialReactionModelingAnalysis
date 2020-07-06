import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns

"""Purpose of script is to open and plot generated histogram files, skipping
the analysis, but instead just putting out the images.

 """

workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData_Stokes\\"
# workingDir = "."
caseName = "TwoInletsTwoColumns_v5.2_ExF"
caseExt = "_histogram\.csv$"
writeMeta = True  # Create a new metadata file
binVel = True  # True to bin velocties, false to skip

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')

for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        data = pd.read_csv(fileName, names=['binID', 'normFreq', 'velVal']