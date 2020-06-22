import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt

""" Purpose of script is to import data for velocity analysis, but will also serve
as a means of doing general data import from comsol. As a reminder, comsol files look like:

%TODO:

-Generate histogram of velocities -> easy?

-Analyze for recirculation zone (how do we tag points as in or not in a zone?)
-Do calculations for residence time?
    How?
-Calcuate the change in pressure over the whole channel
"""

def produceVelPDF(data, nBins):
    """Produce a velocity histogram from the volumes and velocities
    1) Bin the velocities by magnitude
    2) Find the total grid size in each bin
    3) Weight the frequency by multiplying the number of bin entries by the total vol in the bin
    4) Normalize the frequency to make the area under the curve 1. (I think it's sum (freq*area))"""
    # Bin velocities -> is there a better way to do that? We might need it since the bins need to be better at low velocities
    # Option, bin velocities such that the mean and median are within some tolerance or the RSD meets a certain criteria
    # velBin = np.linspace(data.velMag.min(), data.velMag.max(), num=nBins)
    velBin = np.logspace(np.log10(data.velMag.min()), np.log10(data.velMag.max()), num = nBins)
    groups = data.groupby(np.digitize(data.velMag, velBin))
    velVal = groups.velMag.mean()
    # Weight frequencies
    weightedFreq = groups.EleVol.sum()*groups.size()
    # Calculate area under current curve
    totalArea = np.trapz(weightedFreq, velVal)
    normFreq = weightedFreq/totalArea
    return normFreq, velVal, groups, velBin
workingDir = "."
plt.ion()
data = pd.read_table('TestData.txt', header=9, sep='\s+', names=['x', 'y', 'z',
                     'MeshEleNum', 'VolScale', 'EleVol', 'u', 'v', 'w'],
                     usecols=range(0, 9))
avgVol = data['EleVol'].mean()
data['NormScale'] = data['EleVol']/avgVol
data['velMag'] = np.sqrt(data['u']**2+data['v']**2+data['w']**2)
data['velMagScaled'] = data['velMag']/data['NormScale']
data['uScaled'] = data['u']/data['NormScale']
data['vScaled'] = data['v']/data['NormScale']
data['wScaled'] = data['w']/data['NormScale']
# data.hist(column=['velMagScaled', 'uScaled', 'vScaled', 'wScaled'], bins=100)
# data.hist(column=['velMag', 'u', 'v', 'w'], bins=100)
# scat = data.plot.scatter('x', 'y', c='velMag', cmap='hsv', marker='.')
# scat.axis('equal')
normFreq, velVals, velGroups, velBin = produceVelPDF(data, 10000)
plt.figure(4)
plt.plot(velVals, normFreq)
plt.show()
