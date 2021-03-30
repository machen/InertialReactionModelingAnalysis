import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns

""" Purpose of script is to import data for velocity analysis,
but will also serve as a means of doing general data import from comsol.

%TODO:

-Develop a function that can use the extracted parameters to exactly calculate
where the pillar gap is
-Further, you should likely export depth averaged data, and, barring that, depth average
the data you have yourself, likely by coarsening the data.
-Perhaps also output a 2D map of the area considered?


FORMAT IS:
x, y, z, MeshID, MeshVolumeScale, MeshElementVolume, u (m/s), v (m/s), w (m/s),
p, Velocity Magnitude (m/s), Mass flow (kg/s)

DON'T FORGET: X is along channel width, Y is along channel length, Z is along channel depth

-Analyze for recirculation zone (how do we tag points as in or not in a zone?)
-Do calculations for residence time?
    How?

-Metadata file handling:

-Export

Some notes about current results:
The area approximation is kind of inaccurate, as evidenced by studying the math
for a plane just through the main channel width.

For a 0.2 um slice, the cross sectional area is something like 50% of the actual
For a 2 um slice, the area is more like 25% of actual
For a 20 um slice, the area is about 5% of actual

However the total flux is really far off however. I suspect it's because
there is a mismatch beteween the estimated cross sectional area and the
associated flux.

Comparing against the flux in the channel however, the error seems low.
It's easy enough to pull the total flux to see what the rough error is.
"""


def pillarGapCalculation(r1, r2, d):
    """ Use this to properly calculate the true "pillar gap" area depending on
    the available parameters of the model. This calculation will directly assume:
    1) That the edge of the most upstream pillar is at y = 0
    2) That there are only two pillars
    3) Positive y is upstream, negative y is downstream
    4) The simulation base unit is microns

    The gap is defined by the edges of the pillars along the channel length
    and the width of the largest pillar
    """
    if not r2:
        r2 = r1
    xPil = 250
    yPil = -(2*r1+d/2)
    x1 = xPil-max(r1, r2)
    x2 = xPil+max(r1, r2)
    y1 = yPil+d/2
    y2 = yPil-d/2

    return [x1, x2], [y1, y2]


def dataLoader(filename, type='flowdata.txt'):
    if type == 'flowdata.txt':
        data = pd.read_table(fileName, header=9, sep='\s+',
                             names=['x', 'y', 'z', 'meshID', 'eleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow',
                                    'vortX', 'vortY', 'vortZ', 'vortMag'])
    if type == 'chemdata.txt':
        data = pd.read_table(fileName, header=10, sep='\s+',
                             names=['x', 'y', 'z', 'meshID', 'eleVol', 'u',
                                    'v', 'w', 'p', 'velMag', 'massFlow',
                                    'h2o2', 'tcpo', 'cProduct', 'k'])
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


def crossings_nonzero_all(data):
    # Will calculate a zero crossing and return the indices where the crossing occurs
    # Credit to https://stackoverflow.com/questions/3843017/efficiently-detect-sign-changes-in-python
    pos = data > 0
    npos = ~pos
    return ((pos[:-1] & npos[1:]) | (npos[:-1] & pos[1:])).nonzero()[0]


def producePDF(data, nBins=1000, logBin=True, prop="velMag"):
    """
    TO DO: Try to boost computational efficiency by subselecting the data to work on
    (i.e., don't use the columns you don't need)
    subData = data.loc[:,('EleVol', prop)]

    INPUT:
    data: data output from comsol, must include column specified by prop & "eleVol" cols
    nBins: Number of bins you want to use, default is 1000
    logBin: True-> bins are evenly spaced in log space,
    otherwise, linear spacing
    prop: Column name of data you wish to use, default, velMag (velocity magnitude)

    OUTPUT
    normFreq: Normalized frequencies for the defined bins by bin size
    velVal: The mean velocity in each bin
    groups: The binned data group
    velBin: Velocity bin limits used
    Produce a velocity PDF from the input data

    Option, bin velocities such that the mean and median are within
    some tolerance or the RSD meets a certain criteria

    Max velocity is taken as 1.01*max due to the way np.digitize works.
    For a linear bin, without doing this, the last value will not belong to a
    bin that has a proper size definition.

    Alternate bin scheme to enable log binning of negative numbers and ranges
    that cross 0:
    -ID minimum absolute value and maximum absolute value
    -logbin that min to max range
    -duplicate range, reverse it, and then multiply by -1
    -concatenate the two ranges, including 0 between the two
    -pass the ranges directly to the remainder of the calculation
    """

    binMin = data.loc[:, prop].values.min()
    binMax = data.loc[:, prop].values.max()

    if logBin:
        if binMin*binMax < 0:
            # If the binning crosses 0, log binning will not work.
            binPos = data.loc[data.loc[:, prop] > 0, prop].values
            binNeg = data.loc[data.loc[:, prop] < 0, prop].values
            binPosMax = binPos.max()
            binPosMin = binPos.min()
            binNegMax = abs(binNeg).max()
            binNegMin = abs(binNeg).min()
            velBinPos = np.logspace(np.log10(binPosMin*0.99),
                                    np.log10(binPosMax*1.01),
                                    num=int((nBins-2)/2))
            velBinNeg = np.logspace(np.log10(binNegMin*0.99),
                                    np.log10(binNegMax*1.01),
                                    num=int((nBins-2)/2))
            velBinNeg *= -1
            velBinNeg.sort()
            velBin = np.concatenate((velBinNeg, [0], velBinPos))
        elif binMin > 0:
            velBin = np.logspace(np.log10(binMin*0.99),
                                 np.log10(binMax*1.01), num=nBins)
        elif binMin < 0:
            velBin = np.logspace(np.log10(binMax*-.99),
                                 np.log10(binMin*-1.01), num=nBins)
            velBin = np.flip(velBin)
            velBin *= -1
    else:
        if binMin < 0:
            binMin *= 1.01
        else:
            binMin *= 0.99
        if binMax > 0:
            binMax *= 1.01
        else:
            binMax *= 0.99
        if binMin*binMax < 0:
            # If the space crosses 0, insert 0 at the interface
            velBin = np.linspace(binMin, binMax, num=nBins)
            zc = crossings_nonzero_all(velBin)+1
            velBin = np.insert(velBin, zc, [0])

        else:
            velBin = np.linspace(binMin, binMax, num=nBins)
    velVal = (velBin[1:]+velBin[:-1])/2
    normFreq, velBinCopy = np.histogram(data.loc[:, prop].values, bins=velBin,
                                        weights=data.loc[:, 'eleVol'],
                                        density=True)
    return normFreq, velVal, velBin


def extractParams(fileName, nPil=2, caseExt='flowdata.txt'):
    # Produces a dictionary of experimental parameters
    rePat = re.compile('Re(.*?).'+caseExt)
    dPat = re.compile('d(\d+?)_')
    cPat = re.compile('(\d+?)c')
    kPat = re.compile('k(\d+?)_')
    dVal = re.search(dPat, fileName).group(1)
    reVal = re.search(rePat, fileName).group(1)
    if caseExt == 'chemdata.txt':
        cVal = re.search(cPat, fileName).group(1)
        kVal = re.search(kPat, fileName).group(1)
    else:
        cVal = 0
        kVal = 0
    if nPil == 2:
        r1Pat = re.compile('r1_(\d+?)_')
        r2Pat = re.compile('r2_(\d+?)_')
        r1Val = re.search(r1Pat, fileName).group(1)
        r2Val = re.search(r2Pat, fileName).group(1)
        res = {'r1': float(r1Val), 'r2': float(r2Val), 'd': float(dVal),
               'Re': float(reVal), 'c': float(cVal), 'k': float(kVal)}
    if nPil == 1:
        rPat = re.compile('r(\d+?)_')
        rVal = re.search(rPat, fileName).group(1)
        res = {'r1': float(rVal), 'r2': float(rVal), 'd': float(dVal),
               'Re': float(reVal), 'c': float(cVal), 'k': float(kVal)}
    return res


def calcFlowPress(data, params, nu=1.6E-6, c=500E-6, cRatio=0.5,
                  depth=100E-6, regionSize=10):
    """# vel = Re*nu/(c*cRatio)
    # q0 = vel*c*depth
    # nu is default kinematic viscosity for acetonitrile
    # c is the channel width
    # cRatio is the ratio of the inlet channel width to the outlet channel
    # depth is the channel depth
    # Region size is how wide a slice you want to use for pressure averaging
    You should more explicitly calculate pressure of the inflow and outflow
    1) Group data by y value (could be 1 micron intervals) -> data.groupby(data.y)
    Or you could do subDataMin = data.loc[(data.y <= TARGETMAX) & (data.y>= TARGETMIN), :]
    2) find average pressure in each bin which represents a y slice-> WORTH OUTPUTTING POTENTIALLY
    3) Find absolute value of pressure difference of max y bin and min y bin
    4) Output some statistics about that pressure/info about where it's from
    maybe also the distance over which the pressure drop is realized"""
    velInlet = params['Re']*nu/(c*cRatio)
    q0 = velInlet*c*depth  # m^3/s
    yInlet = data.y.values.min()
    inletData = data.loc[(data.y <= yInlet+regionSize), :]
    yOutlet = data.y.values.max()
    outletData = data.loc[(data.y >= yOutlet-regionSize)]
    inletP = inletData.p.values.mean()
    outletP = outletData.p.values.mean()
    deltaP = abs(inletP - outletP)  # Pa
    dx = abs(inletData.y.mean() - outletData.y.mean())
    return deltaP, q0, dx


def calcVortVelAngle(data, uxName, uyName, uzName, wxName, wyName, wzName):
    """each entry is just the string"""
    ux = data.loc[:, uxName].values
    uy = data.loc[:, uyName].values
    uz = data.loc[:, uzName].values
    wx = data.loc[:, wxName].values
    wy = data.loc[:, wyName].values
    wz = data.loc[:, wzName].values
    velMag = np.sqrt(ux**2+uy**2+uz**2)
    vortMag = np.sqrt(wx**2+wy**2+wz**2)
    angle = (ux*wx+uy*wy+uz*wz)/(velMag*vortMag)
    angle = np.degrees(np.arccos(angle))
    return angle


def genOutputFolderAndParams(dataDir, caseName, caseExt, nBins, logBins,
                             binProp, regionName='Result', recircCenter=False,
                             dataRegionX=None, dataRegionY=None,
                             dataRegionZ=None, autoRegion=False):
    if logBins:
        binType = 'log'
    else:
        binType = 'linear'
    outputPath = '..\\{}-{}-{} {} bins\\'.format(regionName,
                                                 binProp, nBins, binType)
    outputFile = outputPath+'Parameters.txt'
    if not os.path.isdir(outputPath):
        os.mkdir(outputPath)
    with open(outputFile, "w") as outFile:
        outFile.write("Data directory: {}\n".format(dataDir))
        outFile.write("Case name: {}\n".format(caseName))
        outFile.write("Case extension: {}\n".format(caseExt))
        outFile.write("Binned property:{}\n".format(binProp))
        outFile.write("Number of bins: {}\n".format(nBins))
        outFile.write("Bin type: {}\n".format(binType))
        outFile.write("X Region: {}\n".format(dataRegionX))
        outFile.write("Y Region: {}\n".format(dataRegionY))
        outFile.write("Z Region: {}\n".format(dataRegionZ))
        outFile.write("Region defined by recirc center: {}\n".format(recircCenter))
        outFile.write("Region defined by geometry: {}\n".format(autoRegion))
    return outputPath


def calcDilutionIndex(data, prop):
    # The values coming out are fine... probably, but we need an appropriate thing to normalize to
    dV = data.eleVol.values
    c = data.loc[:, prop].values
    mTot = np.sum(c*dV)
    mTotMax = np.sum(0.5*dV)
    p = c*dV/mTot
    pMax = 0.5*dV/mTotMax
    e = np.exp(-np.sum(p*np.log(p)))*np.mean(dV)
    vTot = np.sum(dV)
    eMax = np.exp(-np.sum(np.log(pMax)*pMax))*np.mean(dV)
    return e, e/eMax


def recircVolCalc(data, centerCoords, r1, r2, d):
    """
    Calculate the recirculation volume by defining it as a trapezoid contained
    defined by the following points in order:
    -1 and 2 the center of each pillar
    -3 and 4, the intersection of a line passing through the recirculation zone
    center line and the edge of the pillar

      4   3
    1      2
    Trapezoid area therefore is: (abs(x4-x3)+abs(x1-x2))/2*(abs(y_recircCenter-yCenter=250))
    Area of pillar in trapezoid is calculated from the wedge of the circle
    The angle of the trapezoid corner is given by
    cos(theta_2) = abs(x3-x2)/r_2 and cos(theta_1) = abs(x4-x1)/r1
    Subtracted wedge area is then: pi*r_i^2*theta_i/2*pi

    This may actually be unnecessary, since COMSOL doesn't resolve the grid in
    the space that comprises the pillars, so if I define a rectangle that
    runs through the pillar centers and the recirculation center.

    I still need to calculate the recirculation center, which is easy enough to
    figure out if I just sub-select my data to a rough box defined by:

    (xPil1,250), (xPil2, 250), (xPil1,250+r1), (xPil2, 250+r2)
    """
    recircVol = 0
    return recircVol


def estimateFluxes(data, planeWidth=1):
    """# Estimate the mean residence time given:
    data: comsol output, preferrably already pre cut to contain the rough zone of recirculation
    planeWidth: how wide are the planes, which are defined off of coord+/- plane width
    # data should be data that is already pre-selected for the region of interest
    # Break this up into separate functions for calculating:
    center coords
    recircVol -> You should check this against the sum of the eleVol data in the
    recirculation plane, it could easily be that it's just captured

    FIX YOUR RECIRC VOL CALCULATION BY ESTIMATING THE TRAPEZOID AND SUBTRACTING THE
    AREA WHERE THE TRAPEZOID INTERSECTS THE PILLARS
    """
    data = subSelectData(data, xRange=[250, 500], yRange=[-450, -350])  # Use half of the main channel # Should also be selecting for areas that don't intersect the pillar
    midPlane = subSelectData(data, zRange=[50-planeWidth, 50+planeWidth])  # Use the middle plane
    minU = midPlane.velMag.min()
    centerPointRow = midPlane.loc[midPlane.velMag == minU, :]
    centerCoords = [float(centerPointRow.x.values),
                    float(centerPointRow.y.values),
                    float(centerPointRow.z.values)]
    # This is essentially the defined recirculation zone.
    recircData = subSelectData(data, xRange=[250, centerCoords[0]])
    # If we're feeling motivated we should probably ID where flux enters/leaves
    recircVol = recircData.eleVol.sum()  # Native units.
    # This is a cut plane, the total flux through the plane should be 0
    recircPlane = subSelectData(data, xRange=[centerCoords[0]-planeWidth,
                                              centerCoords[0]+planeWidth])
    """Scale the flux off of each elements volume, but assume a constant depth
    sum(velocity*dV/estimated width (constant)) = sum(u*eleVol)/estimatedWidth
    The most accurate would be to determine the actual cross sectional surface
    area of each grid element but that doesn't really seem worth it."""

    # DON'T FORGET XYZ COORDINATES ARE IN um BUT YOU NEED IT IN m!!!!
    # Constant depth
    totalRecircFlux = np.sum(recircPlane.u.values*recircPlane.eleVol.values/(2E-6))
    posRecircPlane = recircPlane.loc[recircPlane.u > 0, :]
    posFlux = np.sum(posRecircPlane.u.values*posRecircPlane.eleVol.values/(2E-6))
    negFlux = totalRecircFlux-posFlux  # Find the opposing flux as well to check accuracy
    return centerCoords, totalRecircFlux, posFlux, negFlux, recircVol

# Read through files in a directory

#workingDir = "..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\RawData\\"
#workingDir = "..\\Comsol5.5\\TwoPillars\\ExF\\ChemData\\RawData"
#workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\FlowData\\RawData\\"
workingDir = "..\\Comsol5.4\\TwoPillars\\Version5\\ExF\\ChemData\\RawData\\"
#workingDir = "TestData"
caseName = "TwoInletsTwoColumns_"
caseExt = "\.chemdata.txt$"
calcFlow = False  # Do Pressure/Flow rate fitting? Only valid with flow
vortAng = False # Calculate the angle between velocity and vorticity vector, will generate data column "angle"
calcChem = True  # Do calculations for PDF from chemistry

print(workingDir)

#PDF Properties

binProp = True  # True to bin velocties, false to skip
dataRegionX = [250, 500]
dataRegionY = [-550, -250]  # [-5000, 250] # Pillar center should be at -400
regionName = 'Pillar Gap Exact'
nBins = 100
logBins = False  # True to use log spaced bins, False to use linear bins
nPil = 2  # Number of pillars in file specification
binProp = 'dCdtNorm'  # Name of column to run PDF on, use 'angle' to do a vort./vel. angle analysis
recircDefinedRegion = False
autoRegion = True

# Chemistry props
diff = 3E-9  # m2/s, H2O2
nu = 1.6E-6  #m^2/s Acetnonitrile kinematic viscosity

# Scipt start

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')

# Check for existence of a metadata file which should be something like:
# caseName+"_meta.csv"
# Purpose of metadata file is to say what we've run already
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2',
                                     'd', 'Re', 'dP', 'q', 'l'])
outFile = genOutputFolderAndParams(workingDir, caseName, caseExt,
                                   nBins, logBins, binProp,
                                   regionName=regionName,
                                   dataRegionX=dataRegionX,
                                   dataRegionY=dataRegionY,
                                   recircCenter=recircDefinedRegion,
                                   autoRegion=autoRegion)
print(outFile)
for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        # Check for fileName already in metaData, skip if so
        data = dataLoader(fileName, type=caseExt[2:-1])
        params = extractParams(fileName, nPil, caseExt=caseExt[2:-1])
        if autoRegion:
            dataRegionX, dataRegionY = pillarGapCalculation(params['r1'], params['r2'], params['d'])
        data = subSelectData(data, xRange=dataRegionX, yRange=dataRegionY)
        params['dP'], params['q'], params['l'] = calcFlowPress(data, params)
        params['recircCenter'], params['totalRecircFlux'], params['posFlux'], params['negFlux'], params['recircVol'] = estimateFluxes(data, 1)
        params['posMRT'] = params['recircVol']/params['posFlux']
        params['negMRT'] = params['recircVol']/abs(params['negFlux'])
        params['fileName'] = fileName
        if vortAng:
            data.loc[:, 'angle'] = calcVortVelAngle(data, 'u', 'v', 'w',
                                                    'vortX', 'vortY', 'vortZ')
        if calcChem:
            kVal = data.loc[data.index[0], 'k']
            dCdt = data.h2o2.values*data.tcpo.values*kVal
            data.loc[:, 'dCdt'] = dCdt
            elementVol = data.eleVol.values
            params['totalVol'] = np.sum(elementVol)
            params['totalProd'] = np.sum(np.multiply(data.cProduct.values, elementVol))
            params['totalTCPO'] = np.sum(np.multiply(data.tcpo.values, elementVol))
            params['totalH2O2'] = np.sum(np.multiply(data.h2o2.values, elementVol))
            params['dCdtAvg'] = np.mean(dCdt)
            params['dCdtStd'] = np.std(dCdt)
            params['dilutionTCPO'], params['reactorTCPO'] = calcDilutionIndex(data, 'tcpo')
            data.loc[:, 'constC'] = data.tcpo.values+data.cProduct.values # Conservative component
            params['dilutionConserv'], params['reactorConserv'] = calcDilutionIndex(data, 'constC')
            dCdtNorm = data.h2o2.values*data.tcpo.values/((params['c']/2)**2)  # Max rate is defined by well mixed, which is going to be cA*cB/(c0/2)**2
            data.loc[:, 'dCdtNorm'] = dCdtNorm
            dCdtMaxNorm = data.h2o2.values*data.tcpo.values/np.max(data.h2o2.values*data.tcpo.values)
            data.loc[:, 'dCdtMaxNorm'] = dCdtMaxNorm

            params['conservative'] = np.sum(data.constC.values*data.eleVol.values)
            params['DH2O2'] = diff
            params['nu'] = nu
            params['velChar'] = params['Re']*nu/500E-6  #Assume 500 um channel width
            params['Pe'] = params['velChar']*params['r1']*2E-6/diff
            params['DaDiff'] = params['k']*params['c']/1000*(2E-6*params['r1'])**2/diff
            params['DaAdv'] = params['k']*params['c']/1000*2E-6*params['r1']/params['velChar']
        if binProp:
            if recircDefinedRegion:
                recircX = [250, params['recircCenter'][0]]
                recircY = [-450+(params['r1']+params['d']/2),
                           -450-(params['r2']+params['d']/2)]
                print('recircX: {}'.format(recircX))
                print('recircY: {}'.format(recircY))
                data = subSelectData(data, xRange=recircX, yRange=recircY)
            normFreq, valMean, valBin = \
                producePDF(data, nBins=nBins, logBin=logBins, prop=binProp)
            pdfData = {'normFreq': normFreq, 'valMean': valMean,
                       'leftBin': valBin[:-1], 'rightBin': valBin[1:]}
            velDF = pd.DataFrame(pdfData)
            velDF.to_csv(outFile+fileName[:-4]+"_histogram.csv")
            plt.figure()
            plt.plot(valMean, normFreq)
            plt.title(binProp)
            plt.xlabel('Average value of bin')
            plt.ylabel('Normalized Frequency (.)')
            plt.savefig(outFile+fileName[:-4]+"_linear.png")
            plt.yscale('log')
            #plt.xscale('log')
            plt.savefig(outFile+fileName[:-4]+"_log.png")
            plt.close()
        metaData = metaData.append(params, ignore_index=True)

metaData.to_csv(outFile+caseName+"_meta.csv")

flowFitData = pd.DataFrame([], columns=['r1', 'r2', 'd', 'linA', 'linB',
                                        'quadA', 'quadB', 'quadC', 'expA',
                                        'expB'])
# Obtain the unique geometries used (i.e. by d, r1, and r2)
uniqueParam = metaData.loc[~metaData.duplicated(['d', 'r1', 'r2']), :]
colorPalette = sns.color_palette('deep', n_colors=len(uniqueParam.index))

f1, ax1 = plt.subplots(1, 1, sharex='col', figsize=(12, 10))
if calcFlow:
    ci = 0
    for i in uniqueParam.index:
        r1 = uniqueParam.loc[i, 'r1']
        r2 = uniqueParam.loc[i, 'r2']
        d = uniqueParam.loc[i, 'd']
        subData = metaData.loc[(metaData.d == d) &
                               (metaData.r1 == r1) & (metaData.r2 == r2), :]
        ax1.plot(subData.q, subData.dP, ls='None', color=colorPalette[ci],
                 marker='o', label='r1 = {}, r2 = {}, d = {}'.format(r1, r2, d))
        linFit = np.polyfit(subData.q, subData.dP, 1)
        quadFit = np.polyfit(subData.q, subData.dP, 2)
        logFit = np.polyfit(np.log(subData.q), np.log(subData.dP), 1)
        interpQ = np.linspace(min(subData.q), max(subData.q), num=100)
        interp_dP = np.polyval(linFit, interpQ)
        interpQuad_dP = np.polyval(quadFit, interpQ)
        interpLog_dP = np.exp(np.polyval(logFit, np.log(interpQ)))
        caseParam = {'r1': r1, 'r2': r2, 'd': d, 'linA': linFit[0],
                     'linB': linFit[1], 'quadA': quadFit[0], 'quadB': quadFit[1],
                     'quadC': quadFit[2], 'expA': logFit[0], 'expB': logFit[1]}
        flowFitData = flowFitData.append(caseParam, ignore_index=True)
        ax1.plot(interpQ, interp_dP, ls='-', color=colorPalette[ci], label="{:.2e}*q+{:2e}".format(*linFit))
        ax1.plot(interpQ, interpQuad_dP, ls='--', color=colorPalette[ci], label="{:.2e}*q^2+{:.2e}*q+{:.2e}".format(*quadFit))
        ax1.plot(interpQ, interpLog_dP, ls='dotted', label="log(dP) = {:.2e}*log(q)+{:.2e}".format(*logFit))
        ci += 1

    flowFitData.to_csv(caseName+"_flowFits.csv")

"""
What do I want to do?

For each unique combination of d, r1, and r2, I want to apply a coloration
and marker to differentiate it, fit Q vs dP, and plot both.

Option 1: Create table of the unique entries
subselect data for each known entry

Option 2: Iterate through d, r1, r2 unique values
Bad because I may not have all combinations.
"""

# Plot deltaP vs Q, also fit a line
# plt.plot(metaData.q, metaData.dP, ls='None', marker='*', color='k',
#          label='Overall')
# linFit = np.polyfit(metaData.q, metaData.dP, 1)
# interp_dP = np.polyval(linFit, metaData.q)
# plt.plot(metaData.q, interp_dP, ls='-', color='k', label='Overall fit')
plt.legend(loc=0)
plt.ylabel('Pressure Difference (Pa)')
plt.xlabel('Flow rate (m^3/s)')
plt.title('Flow rate vs Pressure')
plt.savefig(caseName+'_dPvsQ.png', dpi=600)
