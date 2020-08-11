import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns


def extractParams(fileName):
    # Produces a dictionary of experimental parameters
    r1Pat = re.compile('r1_(\d+?)_')
    r2Pat = re.compile('r2_(\d+?)_')
    rePat = re.compile('Re(.*?).txt')
    dPat = re.compile('d(\d+?)_')
    r1Val = re.search(r1Pat, fileName).group(1)
    r2Val = re.search(r2Pat, fileName).group(1)
    dVal = re.search(dPat, fileName).group(1)
    reVal = re.search(rePat, fileName).group(1)
    res = {'r1': float(r1Val), 'r2': float(r2Val), 'd': float(dVal),
           'Re': float(reVal)}
    return res


fileLoc = "..\\Modeling\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData_FlowOnly\\Pillar region - 500 log bins\\"
fileName = "TwoInletsTwoColumns_v5.1_ExF_FlowOnly_meta.csv"

data = pd.read_csv(fileLoc+fileName, header=0)
uniqueParam = data.loc[~data.duplicated(['d', 'r1', 'r2']), :]

ci = 0
for i in uniqueParam.index:
    r1 = uniqueParam.loc[i, 'r1']
    r2 = uniqueParam.loc[i, 'r2']
    d = uniqueParam.loc[i, 'd']
    subData = data.loc[(data.d == d) &
                           (data.r1 == r1) & (data.r2 == r2), :]
    linFit = np.polyfit(subData.q, subData.dP, 1)
    quadFit = np.polyfit(subData.q, subData.dP, 2)
    logFit = np.polyfit(np.log(subData.q), np.log(subData.dP), 1)
    interpQ = np.linspace(min(subData.q), max(subData.q), num=100)
    interp_dP = np.polyval(linFit, interpQ)
    interpQuad_dP = np.polyval(quadFit, interpQ)
    interpLog_dP = np.exp(np.polyval(logFit, np.log(interpQ)))
    caseParam = {'r1': r1, 'r2': r2, 'd': d, 'linA': linFit[0],
                 'linB': linFit[1], 'quadA': quadFit[0], 'quadB': quadFit[1],
                 'quadC': [2], 'expA': logFit[0], 'expB': logFit[1]}
    ci += 1
