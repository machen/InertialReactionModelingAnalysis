import pandas as pd
import numpy as np
import re

"""Helper functions for use with outputted COMSOL data."""


def dataLoader(fileName, type='flowdata.txt'):
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
    if type == 'velStreamline.txt':
        data = pd.read_table
        data = pd.read_table(fileName, sep='\s+', skiprows=8,
                             names=['x', 'y', 'z', 'sID', 'velMag'])
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