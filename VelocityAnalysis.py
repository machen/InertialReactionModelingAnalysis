import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt

""" Purpose of script is to import data for velocity analysis,
but will also serve as a means of doing general data import from comsol.

%TODO:

-Can I natively capture the column names? -> (\w+? \(.+?\)|\w+?) @ \d+?: will
capture the header name, but it will not capture the resultant parameter info.
Likely that info is not worth it, and should be kept in the file name.

-BUT SHOULD I? No. I should format my data consistently

FORMAT IS:
x, y, z, MeshID, MeshVolumeScale, MeshElementVolume, u (m/s), v (m/s), w (m/s),
p, Velocity Magnitude (m/s), Mass flow (kg/s)

-Generate histogram of velocities -> DONE

-Analyze for recirculation zone (how do we tag points as in or not in a zone?)
-Do calculations for residence time?
    How?
-Calcuate the change in pressure over the whole channel: wait this is easy
it's literally data['p'].max()-data['p'].min(). The flow rate is determined
by the Re number

# vel = Re*nu/(c*c_ratio)
# q0 = vel*c*depth
q0 = Re*nu*depth

-Export: We should output the resulting velocity bins.
"""


def produceVelPDF(data, nBins):
    """Produce a velocity histogram from the volumes and velocities
    1) Bin the velocities by magnitude
    2) Find the total grid size in each bin
    3) Weight the frequency by multiplying the number of bin entries
    by the total vol in the bin
    4) Normalize the frequency to make the area under the curve 1.
    (I think it's sum (freq*area))
    # Bin velocities -> is there a better way to do that?
    We might need it since the bins need to be better at low velocities
    # Option, bin velocities such that the mean and median are within
    some tolerance or the RSD meets a certain criteria
    # velBin = np.linspace(data.velMag.min(), data.velMag.max(), num=nBins)"""
    velBin = np.logspace(np.log10(data.velMag.min()),
                         np.log10(data.velMag.max()), num=nBins)
    groups = data.groupby(np.digitize(data.velMag, velBin))
    velVal = groups.velMag.mean()
    # Weight frequencies
    weightedFreq = groups.EleVol.sum()*groups.size()
    # Calculate area under current curve
    totalArea = np.trapz(weightedFreq, velVal)
    normFreq = weightedFreq/totalArea
    return normFreq, velVal, groups, velBin


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


def calcFlowPress(data, params, nu=1.6E-6, c=500E-6, cRatio=0.5,
                  depth=100E-16):
    # vel = Re*nu/(c*cRatio)
    # q0 = vel*c*depth
    # nu is default kinematic viscosity for acetonitrile
    # c is the channel width
    # cRatio is the ratio of the inlet channel width to the outlet channel
    # depth is the channel depth
    velInlet = params['Re']*nu/(c*cRatio)
    q0 = velInlet*c*depth  # m^3/s
    deltaP = data.p.max()-data.p.min()  # Pa
    return deltaP, q0


workingDir = "."
plt.ion()
nu = 1.6E-6  # m^2/s
fileName = 'TwoInletsTwoColumns_v5.1_Normal_FlowData_d50_r1_50_r2_50_Re0.1.txt'
data = pd.read_table(fileName, header=9, sep='\s+', names=['x', 'y', 'z',
                     'meshID', 'EleVol', 'u', 'v', 'w', 'p', 'velMag',
                     'massFlow'])
avgVol = data['EleVol'].mean()
data['NormScale'] = data['EleVol']/avgVol
data['velMagCalc'] = np.sqrt(np.power(data['u'].values, 2) +
                             np.power(data['v'].values, 2) +
                             np.power(data['w'].values, 2))
data['velMagScaled'] = np.divide(data['velMag'].values,
                                 data['NormScale'].values)
data['uScaled'] = np.divide(data['u'].values, data['NormScale'].values)
data['vScaled'] = np.divide(data['v'].values, data['NormScale'].values)
data['wScaled'] = np.divide(data['w'].values, data['NormScale'].values)
# data.hist(column=['velMagScaled', 'uScaled', 'vScaled', 'wScaled'], bins=100)
# data.hist(column=['velMag', 'u', 'v', 'w'], bins=100)
# scat = data.plot.scatter('x', 'y', c='velMag', cmap='hsv', marker='.')
# scat.axis('equal')
normFreq, velVals, velGroups, velBin = produceVelPDF(data, 1000)
plt.figure(4)
plt.plot(velVals, normFreq)
plt.show()
params = extractParams(fileName)
