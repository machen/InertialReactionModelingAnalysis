from turtle import distance
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns
from itertools import cycle

"""Used for visualizing line plot data exported from comsol.

TODO: Would be nice to setup a module for this project which contains core helper functions
(such as the data loader function, pillar gap calculation, etc.)."""

def pillarGapCalculation(r1, r2, d, includePillar=True, halfRegion=None):
    """ Use this to properly calculate the true "pillar gap" area depending on
    the available parameters of the model.
    This calculation will directly assume:
    1) That the edge of the most upstream pillar is at y = 0
    2) That there are only two pillars
    3) Positive y is upstream, negative y is downstream
    4) The simulation base unit is microns

    The gap is defined by the edges of the pillars along the channel length
    and the width of the largest pillar.

    includePillar either includes or excludes the pillars from the region.

    halfRegion is either 'top', 'bot', or None, indicating to either use
    the top half, bottom half, or to not subdivide the pore throat.
    """
    if not r2:
        r2 = r1
    xPil = 250
    yPil = -(2*r1+d/2)
    if halfRegion == 'top':
        x1 = xPil
        x2 = xPil+max(r1, r2)
    elif halfRegion == 'bot':
        x1 = xPil-max(r1, r2)
        x2 = xPil
    elif not halfRegion:
        x1 = xPil-max(r1, r2)
        x2 = xPil+max(r1, r2)
    if includePillar:
        # Includes the pillar itself
        y1 = yPil+d/2+r1
        y2 = yPil-d/2-r2
    else:
        y1 = yPil+d/2
        y2 = yPil-d/2
    return [x1, x2], [y1, y2]


def dataLoader(filename, type='flowdata.txt'):
    if type == 'flowdata.txt':
        data = pd.read_table(fileName, header=9, sep=r'\s+',
                             names=['x', 'y', 'z', 'meshID', 'eleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow',
                                    'vortX', 'vortY', 'vortZ', 'vortMag'])

    if type == 'chemdata.txt':
        data = pd.read_table(fileName, header=10, sep=r'\s+',
                             names=['x', 'y', 'z', 'meshID', 'eleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow',
                                    'h2o2', 'tcpo', 'cProduct', 'k'])
        # Drop lines where concentration values are less than 0
        data = data.drop(data.loc[(data.h2o2 <= 0)
                         | (data.tcpo <= 0)
                         | (data.cProduct <= 0)].index)
    return data


def subSelectData(data, xRange=None, yRange=None, zRange=None):
    """Assumes that each of the inputs to the function is a tuple containing
     max and min values of x, y, and z that we wish to include. Use for rough
     chopping of the data to exclude main channel flows"""
    if xRange:
        data = data.loc[(data.x > min(xRange)) & (data.x < max(xRange)), :]
    if yRange:
        data = data.loc[(data.y > min(yRange)) & (data.y < max(yRange)), :]
    if zRange:
        data = data.loc[(data.z > min(zRange)) & (data.z < max(zRange)), :]

    return data


def extractParams(fileName, nPil=2, caseExt='flowdata.txt'):
    # Produces a dictionary of experimental parameters
    rePat = re.compile(r'Re(\d+\.?\d*).'+caseExt)
    dPat = re.compile(r'd(\d+\.?\d*)_')
    cPat = re.compile(r'c(\d+\.?\d*)')
    kPat = re.compile(r'k(\d+\.?\d*)_')
    dVal = re.search(dPat, fileName).group(1)
    reVal = re.search(rePat, fileName).group(1)
    if caseExt == 'chemdata.txt':
        cVal = re.search(cPat, fileName).group(1)
        kVal = re.search(kPat, fileName).group(1)
    else:
        cVal = 0
        kVal = 0
    if nPil == 2:
        r1Pat = re.compile(r'r1_(\d+?)_')
        r2Pat = re.compile(r'r2_(\d+?)_')
        r1Val = re.search(r1Pat, fileName).group(1)
        r2Val = re.search(r2Pat, fileName).group(1)
        res = {'r1': float(r1Val), 'r2': float(r2Val), 'd': float(dVal),
               'Re': float(reVal), 'c': float(cVal), 'k': float(kVal)}
    if nPil == 1:
        rPat = re.compile(r'r(\d+?)_')
        rVal = re.search(rPat, fileName).group(1)
        res = {'r1': float(rVal), 'r2': float(rVal), 'd': float(dVal),
               'Re': float(reVal), 'c': float(cVal), 'k': float(kVal)}
    return res


def distanceCalc(data):
    xVal = data.x.values
    yVal = data.y.values
    zVal = data.z.values
    dist = np.sqrt(np.square(xVal[:-1]-xVal[1:])+np.square(yVal[:-1]-yVal[1:])+np.square(zVal[:-1]-zVal[1:]))
    coord = np.cumsum(dist)
    coord = np.insert(coord,0,0)
    data.loc[:,'coord'] = coord
    return data

# SCRIPT INPUTS


workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData-PoreTraverse\\"
caseName = "TwoPillar_v6.*"
caseExt = "\.chemdata.txt$"
nPil = 1
savePlots = False

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
        data = dataLoader(fileName, type=caseExt[2:-1])
        data.sort_values(by=['x','y','z'], inplace=True) # Sort data to make plots reasonable
        data = distanceCalc(data)
        params = extractParams(fileName, nPil, caseExt=caseExt[2:-1])
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