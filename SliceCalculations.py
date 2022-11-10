import datahelper as dh
import os, re
import pandas as pd
import numpy as np

"""Script is meant to do calculate values of things using integrated quantities
(vol. weighted conc., pressure, etc.) in slices, and outputting those results
to a readable spreadsheet

TODO: Implement some kind of warning for setting sliceInput to False. It is better
to just input a file instead of trying to code for all the cases you may run into.

TODO: This doesn't work properly in multiple ways. 1) Using this to
determine an overall volume of a slice results in a mismatch between the dimensions
of the data volume and the sum of the element volumes. 2) Calculated fluid fluxes
at the same positions are incorrect. They do not satisfy mass balance, nor do they
match the integrated fluxes as calculated by COMSOL. DO NOT USE THIS SCRIPT.
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


def sliceVolAvg(dataColumn, eleVol):
    # Use to get a volume weighted average value of something within the slice. Feed a specific column of data values.
    return np.average(dataColumn, weights=eleVol)


def sliceFlux(dataColumn, eleVol, planeWidth=None, conversionFactor=1E-6):
    """Use to calculate the integrated flux of the data col value considering the vol weights.
    If no planeWidth is specified, grid cell i is assumed to be a cube with vol eleVol(i) and equivalent surface area.
    If planeWidth is specified, that is assumed to be the width of a grid cell. conversionFactor should convert the native units
    of the simulation (i.e. x, y, z values) into SI units (m).
    """
    if not planeWidth:
        eleArea = eleVol**(2/3.0)
    else:
        eleArea = eleVol/(planeWidth*conversionFactor)
    flux = dataColumn*eleArea
    return flux.sum()


dataDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\RawData\\"
workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\SliceData\\"
caseName = 'TwoPillar_v6_ExF_c3_k2000_r100_d100_Re100'
caseExt = 'chemdata.txt'
nPil = 1  # Should be the number of pillars specified in the file name
sliceInput = "TestSlices.csv"  # Set true to use a parameter file containing names and slices to use for a given geometry. Otherwise, script should interpolate based on file geometry
planeWidth = 5  # Set to larger if you have empty planes, in the native unit of the simulation.
useWidth = True  # Set to false if you want to estimate elements as cubes

filePat = re.compile(caseName+'.*?'+caseExt+'$')
if not os.path.exists(workingDir):
    os.mkdir(workingDir)
dataList = os.listdir(dataDir)
if sliceInput:
    # load slice list from file. Useful for cases with the same geometry
    try:
        slices = pd.read_csv(workingDir+sliceInput, header=0)
    except FileNotFoundError:
        print('No slice specification, stopping script.')
        raise
else:
    sliceTemp = pd.read_csv('SliceTemplate.csv', header=0)


# Iterate through files, picking the ones matching the caseName and extension
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2', 'd', 'Re',
                                     'kVal', 'c', 'sliceName', 'range1Min',
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
            out['planeWidth'] = planeWidth
            out['useWidth'] = useWidth
            # Switch for estimating element volume depth
            if useWidth:
                sliceWidth = planeWidth
            else:
                sliceWidth = None
            dataSlice = dh.slicePlane(data,
                                      [slice['range1Min'], slice['range1Max']],
                                      [slice['range2Min'], slice['range2Max']],
                                      slice['loc'], slice['type'], planeWidth)
            eleVol = dataSlice.eleVol
            if slice['type'] == 'xy':
                normalVel = dataSlice.w
            elif slice['type'] == 'yz':
                normalVel = dataSlice.u
            elif slice['type'] == 'xz':
                normalVel = dataSlice.v
            # Calc volume weighted pressure.
            out['volAvgPress'] = sliceVolAvg(dataSlice.p, eleVol)
            # Calc volumetric fluxes, fluid, solute masses
            out['fluidFlux'] = sliceFlux(normalVel, eleVol, sliceWidth)
            if caseExt == 'chemdata.txt':
                out['tracerTCPOFlux'] = sliceFlux((dataSlice.tcpo +
                                                   dataSlice.cProduct) *
                                                  normalVel,
                                                  eleVol, sliceWidth)
                out['tcpoFlux'] = sliceFlux(dataSlice.tcpo *
                                            normalVel,
                                            eleVol, sliceWidth)
                out['h2o2Flux'] = sliceFlux(dataSlice.h2o2 *
                                            normalVel, eleVol, sliceWidth)
                out['prodFlux'] = sliceFlux(dataSlice.cProduct *
                                            normalVel, eleVol,
                                            sliceWidth)
            out['massFlux'] = dataSlice.massFlow.sum()
            out['sliceVol'] = eleVol.sum()
            metaData = metaData.append(out, ignore_index=True)

if sliceInput:
    metaData.to_csv(workingDir+sliceInput[:-4]+"_processed.csv")
else:
    metaData.to_csv(workingDir+'AutoSlice.csv')
