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
data.hist(column=['velMagScaled', 'uScaled', 'vScaled', 'wScaled'], bins=100)
data.hist(column=['velMag', 'u', 'v', 'w'], bins=100)
plt.show()
