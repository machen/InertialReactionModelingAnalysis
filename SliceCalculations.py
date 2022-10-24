import datahelper as dh
import os, re
import pandas as pd
import numpy as np

"""Script is meant to do calculate values of things using integrated quantities
(vol. weighted conc., pressure, etc.) in slices, and outputting those results
to a readable spreadsheet

TODO: Implement some kind of warning for setting sliceInput to False. It is better
to just input a file instead of trying to code for all the cases you may run into.

TODO: Implement sliceInput as a file instead of bool
"""


def genSliceFromPillar(data, xPil, yPil, rPil, half=None):
    # Should generate coordinates that can be used to get slices for an XZ plane
    # Slice can be around the pillar, or inlcude the
    zMax = data.z.max()
    zMin = data.z.min()
    xMax = data.x.max()
    xMin = data.x.min()
    if not half:
        # Assumes you want the whole sim just sliced at one line
        slice = [xMin, xMax, zMin, zMax, yPil, 'xz']
        return slice
    elif half == 'top':
        # Slices above and below the pillar, only works for a single line of pillars
        topSlice = [xPil+rPil, xMax, zMin, zMax, yPil, 'xz']
        return topSlice
    elif half == 'bot':
        botSlice = [xMin, xPil-rPil, zMin, zMax, yPil, 'xz']
        return botSlice
    else:
        return


def gen2PillarSlices(data, slices, r1, r2, d):
    # Assumes basic 2 pillar geometry which is used by dh.PillarGapCalculation
    # Creates slices at the top and bottom of the pillars as well as start
    # and end of the pore throat
    xPil = 250  # Assumes symmetric pillars
    yPil1 = 0-r1  # Asumes 1st pillar edge is at y=0
    yPil2 = -(2*r1+d+r2)
    slices.loc[0, 'sliceName'] = 'p1Top'
    slices.iloc[0, 1:] = genSliceFromPillar(data, xPil, yPil1, r1, half='top')
    slices.loc[1, 'sliceName'] = 'p1Bot'
    slices.iloc[1, 1:] = genSliceFromPillar(data, xPil, yPil1, r1, half='bot')
    slices.loc[2, 'sliceName'] = 'p2Top'
    slices.iloc[2, 1:] = genSliceFromPillar(data, xPil, yPil2, r2, half='top')
    slices.loc[3, 'sliceName'] = 'p2Bot'
    slices.iloc[3, 1:] = genSliceFromPillar(data, xPil, yPil2, r2, half='bot')
    slices.loc[4, 'sliceName'] = 'throatStart'
    slices.iloc[4, 1:] = genSliceFromPillar(data, xPil, 0, r1)
    slices.loc[5, 'sliceName'] = 'throatEnd'
    slices.iloc[5, 1:] = genSliceFromPillar(data, xPil, yPil2-r2, r2)
    return slices


dataDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\RawData\\"
workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\SliceData\\"
caseName = 'TwoPillar_v6'
caseExt = 'chemdata.txt'
nPil = 1  # Should be the number of pillars specified in the file name
# TODO: Make sliceInput a filename
sliceInput = False  # Set true to use a parameter file containing names and slices to use for a given geometry. Otherwise, script should interpolate based on file geometry


filePat = re.compile(caseName+'.*?'+caseExt+'$')
if not os.path.exists(workingDir):
    os.mkdir(workingDir)
dataList = os.listdir(dataDir)
workingList = os.listdir(workingDir)
if sliceInput:
    # load slice list from file. Useful for cases with the same geometry
    try:
        sliceTemp = pd.read_csv(workingDir+'Slices.csv', header=0)
    except FileNotFoundError:
        print('No slice specification, stopping script.')
        raise
else:
    sliceTemp = pd.read_csv('SliceTemplate.csv', header=0)


# Iterate through files, picking the ones matching the caseName and extension
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2', 'd', 'Re',
                                     'k', 'c', 'sliceName', 'range1Min',
                                     'range1Max', 'range2Min', 'range2Max',
                                     'loc', 'type'])
for fileName in dataList:
    if re.match(filePat, fileName):
        print(fileName)
        data = dh.dataLoader(dataDir+fileName, caseExt)
        params = dh.extractParams(fileName, nPil, caseExt)
        params['fileName'] = fileName  # TODO: This should be built into extractParams()
        if not sliceInput:
            # Create slice list from geometry of pillars. Need helper function to calculate pillar locations from parameters.
            if nPil == 1:
                slices = gen2PillarSlices(data, sliceTemp, params['r1'],
                                          params['r2'], params['d'])
            else:
                raise ValueError('Implied file geometry is not implemented.')
        # Defines the data representing an individual file and slice
        out = params
        for indexVal in slices.index:
            slice = slices.loc[indexVal, :]
            # Put slice info into params for export
            out.update(slice.to_dict())
            dataSlice = dh.slicePlane(data,
                                      [slice['range1Min'], slice['range1Max']],
                                      [slice['range2Min'], slice['range2Max']],
                                      slice['loc'], slice['type'])
            # Calc volume weighted pressure.
            # Calculate ALL vol. weighted params
            wDataSliceSum = dataSlice.mul(dataSlice.eleVol, 0).sum()
            wDataSliceSum.index += 'VolWeight'
            dataSliceSum = dataSlice.sum()
            # Included summed parameters, means we don't need dataformat
            out.update(dataSliceSum.to_dict())
            out.update(wDataSliceSum.to_dict())
            metaData = metaData.append(out, ignore_index=True)

if sliceInput:
    metaData.to_csv(workingDir+'processedSlices.csv')
else:
    metaData.to_csv(workingDir+'AutoSlice.csv')
