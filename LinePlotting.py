import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns
from itertools import cycle
import datahelper as dh

"""Used for visualizing line plot data exported from comsol.
"""


def distanceCalc(data):
    xVal = data.x.values
    yVal = data.y.values
    zVal = data.z.values
    dist = np.sqrt(np.square(xVal[:-1]-xVal[1:])+np.square(yVal[:-1]-yVal[1:])+np.square(zVal[:-1]-zVal[1:]))
    coord = np.cumsum(dist)
    coord = np.insert(coord, 0, 0)
    data.loc[:,'coord'] = coord
    return data

# SCRIPT INPUTS


workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData-PoreTraverse\\"
caseName = "TwoPillar_v6.*"
caseExt = "\.chemdata.txt$"
nPil = 1
savePlots = False
diffPlot = False

# Script here

sns.set_context('poster', font_scale=1)
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Cambria'

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')

f1, ax1 = plt.subplots(1,1,sharex='col',figsize=(12,10)) #H2O2
f2, ax2 = plt.subplots(1,1,sharex='col',figsize=(12,10)) #TCPO
f3, ax3 = plt.subplots(1,1,sharex='col',figsize=(12,10)) #PRODUCT
f4, ax4 = plt.subplots(1,1,sharex='col',figsize=(12,10)) #TRACER

lsCycle = cycle(['-', '--', '-.',':'])
colorCycle = cycle(['b','tab:orange','g','r','tab:pink','k'])
dVals = {}
reVals = {}


for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        data = dh.dataLoader(fileName, type=caseExt[2:-1])
        data.sort_values(by=['x','y','z'], inplace=True) # Sort data to make plots reasonable
        data = distanceCalc(data)
        params = dh.extractParams(fileName, nPil, caseExt=caseExt[2:-1])
        if params['d'] not in dVals:
            dVals[params['d']] = next(lsCycle)
        ls = dVals[params['d']]
        if params['Re'] not in reVals:
            reVals[params['Re']] = next(colorCycle)
        color = reVals[params['Re']]
        ax1.plot(data.coord, data.h2o2, ls=ls, color=color,
                 label='H2O2, ReChan: {Re}, d: {d}'.format(Re=params['Re'],
                                                           d=params['d']))
        ax2.plot(data.coord, data.tcpo, ls=ls, color=color,
                 label='TCPO, ReChan: {Re}, d: {d}'.format(Re=params['Re'],
                                                           d=params['d']))
        ax3.plot(data.coord, data.cProduct, ls=ls, color=color,
                 label='H2O2, ReChan: {Re}, d: {d}'.format(Re=params['Re'],
                                                           d=params['d']))
        ax4.plot(data.coord, data.tcpo+data.cProduct, ls=ls, color=color,
                 label='Tracer, ReChan: {Re}, d: {d}'.format(Re=params['Re'],
                                                             d=params['d']))
        # Plot all trends for a given Re
        data.plot(x='coord', y=['h2o2', 'tcpo', 'cProduct'], legend=True,
                  xlabel='Coordinate', ylabel='Concentration', figsize=(12,10))
        plt.title('ReChan: {Re}, d: {d}'.format(Re=params['Re'],
                                                d=params['d']))
        sns.despine()
        if savePlots:
            plt.savefig(fileName[:-3]+'Conc.svg', dpi=600, format='svg')
        if diffPlot:
            diff = data.diff()
            diff.loc[0, :] = 0
            data.loc[:, 'dtcpo'] = diff.tcpo/diff.coord
            data.loc[:, 'dh2o2'] = diff.h2o2/diff.coord
            data.loc[:, 'dprod'] = diff.cProduct/diff.coord
            data.plot(x='coord', y=['dh2o2','dtcpo','dprod'], legend=True,
                    xlabel='Coordinate', ylabel=('dConc/dCoord'), figsize=(12,10))
            plt.title('ReChan: {Re}, d: {d}'.format(Re=params['Re'],
                                                    d=params['d']))
            sns.despine()
            if savePlots:
                plt.savefig(fileName[:-3]+'Conc.svg', dpi=600, format='svg')
            plt.savefig(fileName[:-3]+'dConc.svg', dpi=600, format='svg')


ax1.set_xlabel('Coordinate')
ax1.set_ylabel('Concentration')
ax1.set_title('H2O2')
ax1.legend(loc=0)
sns.despine(f1)


ax2.set_xlabel('Coordinate')
ax2.set_ylabel('Concentration')
ax2.set_title('TCPO')
ax2.legend(loc=0)
sns.despine(f2)


ax3.set_xlabel('Coordinate')
ax3.set_ylabel('Concentration')
ax3.set_title('Product')
ax3.legend(loc=0)
sns.despine(f3)


ax4.set_xlabel('Coordinate')
ax4.set_ylabel('Concentration')
ax4.set_title('Tracer')
ax4.legend(loc=0)
sns.despine(f4)

if savePlots:
    f1.savefig('H2O2.svg', dpi=600, format='svg')
    f2.savefig('TCPO.svg', dpi=600, format='svg')
    f3.savefig('Product.svg', dpi=600, format='svg')
    f4.savefig('TCPO Tracer.svg', dpi=600, format='svg')

plt.ion()
plt.show()