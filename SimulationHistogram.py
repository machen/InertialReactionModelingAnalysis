import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import os
import seaborn as sns
import datahelper as dh
# from mpl_toolkits.mplot3d import Axes3D


""" Purpose of script is to import data for velocity analysis,
but will also serve as a means of doing general data import from comsol.

%TODO:

-Further, you should likely export depth averaged data, and, barring that,
 depth average the data you have yourself, likely by coarsening the data.
-Perhaps also output a 2D map of the area considered?

Is there a way to boost the speed of the readin?


FORMAT IS:
x, y, z, MeshID, MeshVolumeScale, MeshElementVolume, u (m/s), v (m/s), w (m/s),
p, Velocity Magnitude (m/s), Mass flow (kg/s)

DON'T FORGET: X is along channel width, Y is along channel length, Z is along
channel depth

-Analyze for recirculation zone (how do we tag points as in or not in a zone?)
-Do calculations for residence time?
    How?

-Metadata file handling:

-Export

Some notes about current results:
The area approximation is kind of inaccurate, as evidenced by studying the math
for a plane just through the main channel width.

For a 0.2 um slice, the cross sectional area is ~50% of the actual
For a 2 um slice, the area is more like 25% of actual
For a 20 um slice, the area is about 5% of actual

However the total flux is really far off however. I suspect it's because
there is a mismatch beteween the estimated cross sectional area and the
associated flux.

Comparing against the flux in the channel however, the error seems low.
It's easy enough to pull the total flux to see what the rough error is.
"""


def crossings_nonzero_all(data):
    """Will calculate a zero crossing and
    return the indices where the crossing occurs
    Credit to https://stackoverflow.com/questions/3843017/
    efficiently-detect-sign-changes-in-python"""
    pos = data > 0
    npos = ~pos
    return ((pos[:-1] & npos[1:]) | (npos[:-1] & pos[1:])).nonzero()[0]


def producePDF(data, nBins=1000, logBin=True, prop="velMag"):
    """
    TODO: Try to boost computational efficiency by subselecting the data
    to work on
    (i.e., don't use the columns you don't need)
    subData = data.loc[:,('EleVol', prop)]

    INPUT:
    data: data output from comsol,
    must include column specified by prop & "eleVol" cols
    nBins: Number of bins you want to use, default is 1000
    logBin: True-> bins are evenly spaced in log space,
    otherwise, linear spacing
    prop: Column name of data you wish to use, default,
    velMag (velocity magnitude)

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
    1) Group data by y value (could be 1 micron intervals)
    -> data.groupby(data.y)
    Or you could do subDataMin =
    data.loc[(data.y <= TARGETMAX) & (data.y>= TARGETMIN), :]
    2) find average pressure in each bin which represents a y slice->
    WORTH OUTPUTTING POTENTIALLY
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
                             dataRegionZ=None, autoRegion=False,
                             includePillar=False,
                             halfRegion=None):
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
        outFile.write("Region defined by recirc center: {}\n"
                      .format(recircCenter))
        outFile.write("Region defined by geometry: {}\n".format(autoRegion))
        if halfRegion:
            outFile.write("Using {} half\n".format(halfRegion))
        outFile.write("Region includes pillars: {}\n".format(includePillar))
    return outputPath


def calcDilutionIndex(data, prop):
    """ Calculates the dilution index following Kitidanis, 1994.
    Note that we use the discrete case, not the continuous.
    Notably, the continuous case has a thorny issue where the units
    don't match that of the discrete case."""
    dV = data.eleVol.values
    c = data.loc[:, prop].values
    mTot = np.sum(c*dV)
    """ If we define well mixed as the concentration is the same everywhere,
    the contribution of the "c" term will normalize out
    """
    mTotMax = np.sum(dV)
    p = c*dV/mTot
    # Same as above
    pMax = dV/mTotMax
    e = np.exp(-np.sum(p*np.log(p)))*np.mean(dV)
    eMax = np.exp(-np.sum(np.log(pMax)*pMax))*np.mean(dV)
    return e, e/eMax


def centerPointEstimation(data, planeWidth=1,
                          r1=100, r2=100, d=100, useMid=True):
    """Given a dataset, should give the estimated centerpoint using the
    point with the lowest velocity in the middle (z=50) plane of the channel.
    """
    # Calculate the exact space between the pillars
    xGap, yGap = dh.pillarGapCalculation(r1, r2, d, False)
    xGap = [260, xGap[1]-10]  # Back off of edges by 10 um
    yGap = [yGap[0]-10, yGap[1]+10]  # Back off of edges by 10 um
    data = dh.subSelectData(data, xRange=xGap, yRange=yGap, zRange=[10, 90])
    if useMid:
        midPlane = dh.subSelectData(data,
                                 zRange=[50-planeWidth*0.5, 50+planeWidth*0.5])
        if midPlane.empty:
            planeWidthAdjust = planeWidth
            while midPlane.empty:
                planeWidthAdjust += planeWidth
                print("Plane width too small, incrementing planeWidth")
                midPlane = dh.subSelectData(data,
                                         zRange=[50-planeWidthAdjust*0.5,
                                                 50+planeWidthAdjust*0.5])
                if planeWidthAdjust >= 10*planeWidth:
                    print('WARNING: Planewidth too large, writing null')
                    return None
        minU = midPlane.velMag.min()
        centerPointRow = midPlane.loc[midPlane.velMag == minU, :]
    else:
        minU = data.velMag.min()
        centerPointRow = data.loc[data.velMag == minU, :]
    centerCoords = [float(centerPointRow.x.values),
                    float(centerPointRow.y.values),
                    float(centerPointRow.z.values)]
    return centerCoords


def estimateRecircFlux(recircData, centerCoords, recircVol, planeWidth=1):
    """Given a set of data that contains a SINGLE recirculation zone
    and the center coordinates for that zone, calculates the recirculation
    flux through the zone by estimating the flux that cuts through the
    centerplane """
    recircPlane = dh.subSelectData(recircData, xRange=[
                                centerCoords[0]-planeWidth*0.5,
                                centerCoords[0]+planeWidth*0.5])
    totalRecircFlux = np.sum(recircPlane.u.values*recircPlane.eleVol.values/(planeWidth*1E-6))
    posRecircPlane = recircPlane.loc[recircPlane.u > 0, :]
    posFlux = np.sum(posRecircPlane.u.values*posRecircPlane.eleVol.values/(planeWidth*1E-6))
    negFlux = totalRecircFlux-posFlux  # Find the opposing flux as well to check accuracy
    return totalRecircFlux, posFlux, negFlux


def selectRecircZoneBasic(data, r1, r2, d, includePillar):
    """ Does a basic selection of recirculation zone, defining it by:
    The pillar edges, the channel centerline, and the farthest point where a negative
    x velocity still exists.
    """
    xEdge = max(data.loc[data.loc[:, 'v'] > 0, 'x'])
    recircData = dh.subSelectData(data, xRange=[250, xEdge])
    xGap, yGap = dh.pillarGapCalculation(r1, r2, d, includePillar)
    recircData = dh.subSelectData(recircData, yRange=yGap)
    return recircData, xEdge


def selectRecircZoneNegVel(data, r1, r2, d, gridSize=1):
    """ As opposed to subSelectData which is just a basic square cut, this
    is supposed to very accurately delineate the recirculation zone by defining
    the following boundaries:
    1) The middle line of the device aligned with the flow axis (x=250)
    2) The middle of the left pillar
    3) The middle of the right pillar
    4) A curve delineated by the points which have negative flow velocity that
    are farthest from the line in 1)
    4 is the most difficult part of the algorithm, since this is an actual curve
    that depends on the actual system.digitize

    Algorithm:
    Divide the yRange into spots defined by grid resolution (numpy)
    Find the point farthest from the centerline which still has negative velocity in  # noqa: E501
    each subspot
    Select or drop data outside of that range

    The question is can I do this in a more...optimal way.

    For example, could I pre bin the data based on the grid resolution and
    then just get the min from that?
    for example:
    bins = np.linspace(data.y.min(),data.y.max(),num=gridSize)
    data['Ybin'] = pd.cut(data.y,bins,label=bins[:-1])
    grouped = data.groupby(by='Ybin')
    xVals = grouped.v.min()
    data['xLim'] = xVals # This doesn't map right, but this is the idea of what I want
    # Criterion for retention is that y is in yrange and x < xVal[yRange]
    data.drop(data.x >data.xVals)

    """
    # Boundaries 2 and 3 covered by this
    xRange, yRange = dh.pillarGapCalculation(r1, r2, d, True)
    data = dh.subSelectData(data, xRange=xRange, yRange=yRange)
    # Covers middle line in 1)
    data = dh.subSelectData(data, xRange=[250, 500])
    bins = np.arange(data.y.min(), data.y.max(), step=gridSize)
    # bin data, label is left end
    data['yBin'] = pd.cut(data.y, bins, labels=bins[:-1])
    for leftBin in bins[:-1]:
        # Subselect the "grid" block
        subData = data.loc[data.yBin == leftBin, :]
        if not (subData.v>0).any():
            data = data.loc[data.yBin!=leftBin, :]  # If data has no opposing velocities, remove and move on
            continue
        # Apparently upgrading to 1.2.4 pandas broke this. FUN FUN FUN
        # Seems like that there are "missing labels" in the subdataarg1, arg2
        maxX = max(subData.loc[subData.v > 0, 'x'])
        data.drop(data.loc[(data.yBin == leftBin) & (data.x>maxX),:].index,
        inplace=True) # Slice out data in the grid bock
    return data


def calcRecircBoundsFluxInt(data, gridSize=1):
    """ Strictly speaking we should be integrating in x from 250 to some point
    rather than integrating from the top and the bottom because we "know" that the
    separating surface should be roughly in the yz axis

    This requires that you submit a data set which contains only a single
    recirculation zone (i.e., running from 250 to 500 or 0 to 250)
    """
    recircData = data.copy()
    # Define grid edges to rebin the data into
    gridX = np.arange(round(recircData.x.min()), round(recircData.x.max())+gridSize, step=gridSize)
    gridY = np.arange(round(recircData.y.min()), round(recircData.y.max())+gridSize, step=gridSize)
    gridZ = np.arange(round(recircData.z.min()), round(recircData.z.max())+gridSize, step=gridSize)
    # ID Grid elements by the centerpoint of the bin
    gridXCoords = (gridX[1:]+gridX[:-1])/2
    gridYCoords = (gridY[1:]+gridY[:-1])/2
    gridZCoords = (gridZ[1:]+gridZ[:-1])/2
    recircData['xCoord'] = pd.cut(recircData.x, gridX, right=False, labels=gridXCoords)
    recircData['yCoord'] = pd.cut(recircData.y, gridY, right=False, labels=gridYCoords)
    recircData['zCoord'] = pd.cut(recircData.z, gridZ, right=False, labels=gridZCoords)
    boundCoords = np.zeros((len(gridYCoords)*len(gridZCoords),3))
    index = 0
    # Raster over each grid block. Would nice to vectorize this.
    for yCoord in gridYCoords:
        for zCoord in gridZCoords:
            # Pick the column of data corresponding to the y and z coords we're picking
            subData = recircData.loc[(recircData.yCoord == yCoord) & (recircData.zCoord == zCoord), :].copy()
            # Sort the data from low to high (i.e. from 250 up, may need to smartly detect where it changes)
            subData.sort_values(by='xCoord', ascending=True, inplace=True)
            xCoords = subData.xCoord.values
            u = subData.v.values # Should be v because y direction is along channel
            vol = subData.eleVol.values
            flux = u*np.square(np.cbrt(vol))
            cumFlux = np.cumsum(flux, axis=0)
            indices = crossings_nonzero_all(cumFlux)
            if len(indices) == 0:
                recircData = recircData.drop(recircData.loc[(recircData.yCoord == yCoord) & (recircData.zCoord == zCoord),:].index)
                boundCoords[index, :] = [0, yCoord, zCoord]
                index +=1
                continue
            xCoord = xCoords[min(indices)] # This is absolutely hacky and you should make something more rigorous
            recircData = recircData.drop(recircData.loc[(recircData.yCoord == yCoord) & (recircData.zCoord == zCoord) & (recircData.xCoord > xCoord),:].index)
            boundCoords[index, :] = [xCoord, yCoord, zCoord]
            index += 1
    return recircData, boundCoords

def estimateFluxes(data, planeWidth=1, r1=100, r2=100, d=100, useMid=True):
    """
    DEFINITELY BREAK UP INTO SEPARATE FUNCTIONS.
    # Estimate the mean residence time given:
    data: comsol output, preferrably already pre cut to contain the rough zone of recirculation
    planeWidth: how wide are the planes, which are defined off of coord+/- plane 0.5*width
    r1: the first pillar radius in microns
    r2: the second pillar radius in microns
    r3: the gap distance, in microns
    # data should be data that is already pre-selected for the region of interest
    # Break this up into separate functions for calculating:
    center coords
    recircVol -> You should check this against the sum of the eleVol data in the
    recirculation plane, it could easily be that it's just captured

    An issue with this is that it will too easily try to capture values right near
    the midline, or right against the pillars. The easy solution is to just... back off of those a little.

    We can use the pillarGapCalculation to generate the exact space between the pillars
    """
    data = dh.subSelectData(data, xRange=[250, 500])  # Use half of the main channel # Should also be selecting for areas that don't intersect the pillar
    xGap, yGap = pillarGapCalculation(r1, r2, d)  # Calculate the exact space between the pillars
    xGap = [260, xGap[1]-10]  # Back off of edges by 10 um
    yGap = [yGap[0]-10, yGap[1]+10]  # Back off of edges by 10 um
    if useMid:
        midPlane = dh.subSelectData(data, xRange=xGap, yRange=yGap,
                                 zRange=[50-planeWidth*0.5, 50+planeWidth*0.5]) # Use the middle plane
        if midPlane.empty:
            planeWidthAdjust = planeWidth
            while midPlane.empty:
                planeWidthAdjust += planeWidth
                print("Plane width too small, incrementing planeWidth")
                midPlane = dh.subSelectData(data, xRange=xGap, yRange=yGap,
                                         zRange=[50-planeWidthAdjust*0.5, 50+planeWidthAdjust*0.5])
                if planeWidthAdjust >= 10*planeWidth:
                    print('WARNING: Planewidth too large, writing null')
                    return [0, 0, 0], 0, -1, -1, 0, 0
        minU = midPlane.velMag.min()
    else:
        midU = data.velMag.min()
    centerPointRow = midPlane.loc[midPlane.velMag == minU, :]
    centerCoords = [float(centerPointRow.x.values),
                    float(centerPointRow.y.values),
                    float(centerPointRow.z.values)]
    xEdge = max(data.loc[data.loc[:, 'v'] > 0, 'x'])  # serach for farthest point where flow goes against mean direction
    # This is essentially the defined recirculation zone.
    recircData = dh.subSelectData(data, xRange=[250, xEdge])
    # If we're feeling motivated we should probably ID where flux enters/leaves
    recircVol = recircData.eleVol.sum()  # Native units.
    # This is a cut plane, the total flux through the plane should be 0
    recircPlane = dh.subSelectData(data, xRange=[centerCoords[0]-planeWidth*0.5,
                                              centerCoords[0]+planeWidth*0.5])
    """Scale the flux off of each elements volume, but assume a constant depth
    sum(velocity*dV/estimated width (constant)) = sum(u*eleVol)/estimatedWidth
    The most accurate would be to determine the actual cross sectional surface
    area of each grid element but that doesn't really seem worth it."""

    # Assumes the plane is uniform width, thus scaling eleVol by the width will give the surface area for flux
    totalRecircFlux = np.sum(recircPlane.u.values*recircPlane.eleVol.values/(planeWidth*1E-6))
    posRecircPlane = recircPlane.loc[recircPlane.u > 0, :]
    posFlux = np.sum(posRecircPlane.u.values*posRecircPlane.eleVol.values/(planeWidth*1E-6))
    negFlux = totalRecircFlux-posFlux  # Find the opposing flux as well to check accuracy
    return centerCoords, totalRecircFlux, posFlux, negFlux, recircVol, xEdge


def plotDataSet(data, label):
    # Plot the selected data.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(data.x, data.y, data.z,c=data.x, label=label, marker='.')
    ax.set_xlabel('x (Cross flow direction)')
    ax.set_ylabel('y (Mean flow direction)')
    ax.set_zlabel('z')
    ax.legend()
    return ax


def plotPoints(ax, coords):
    ax.scatter(coords[0], coords[1], coords[2], color='r')
    return

"""SCRIPT INPUTS"""

#workingDir = "..\\Comsol5.5\\TwoPillars\\ExF\\FlowDatawVorticity\\RawData\\"
# workingDir = "..\\Comsol5.5\\TwoPillars\\ExF\\ChemData\\RawData"
workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\ChemData\\RawData\\"
# workingDir = "..\\Comsol5.4\\TwoPillars\\Version6\\ExF\\FlowData\\RawData\\"
# workingDir = "TestData"
# workingDir = "..\\Working Data\\ChemData\\RawData\\"
caseName = "TwoPillar_v6"
# caseName = "TwoPillar_v6_ExF_FlowOnly_r100_d100_Re100"
caseExt = r"\.chemdata.txt$"
# caseExt = "\.flowdata.txt$"
calcFlow = False  # Do Pressure/Flow rate fitting? Only valid with flow
vortAng = False  # Calculate the angle between velocity and vorticity vector, will generate data column "angle"
calcChem = True  # Do calculations for PDF from chemistry

print(workingDir)

#PDF Properties

testMode = False  # Set to true to use only one file, which you have to specify
plotData = False

binProp = True  # True to bin values defined by binProp, false to skip
dataRegionX = [150, 350]
dataRegionY = [-5000, 250]
useMid = True  # Use middle plane for calculating recirc center?
regionName = 'Pillar gap pillar exclusive'
nBins = 100
logBins = False  # True to use log spaced bins, False to use linear bins
nPil = 1  # Number of pillars in file specification
binProp = 'dCdt'  # Name of column to run PDF on, use 'angle' to do a vort./vel. angle analysis
estimateRecircCenter = False
recircDefinedRegion = False  # Will cut data to strictly defined single recirculation zone (x=250+)
autoRegion = True
halfRegion = None # 'top', 'bot', or None (uses whole region)
includePillar = False
maxValue = 4.39  #  4.39 for dC/dt sim, 100 um pillar gap. 3 for TCPO/product sims. User input value for calculating dCdtMaxNorm, this should be drawn from the highest observed value in simulated cases
metaData = pd.DataFrame([], columns=['fileName', 'r1', 'r2',
                                     'd', 'Re', 'dP', 'q', 'l'])
# Chemistry props
diff = 3E-9  # m2/s, H2O2
nu = 4.3E-7  #m^2/s Acetnonitrile kinematic viscosity


""" Script start """

os.chdir(workingDir)
filePat = re.compile(caseName+'.*?'+caseExt)
fileList = os.listdir('.')

# Check for existence of a metadata file which should be something like:
# caseName+"_meta.csv"
# Purpose of metadata file is to say what we've run already

outFile = genOutputFolderAndParams(workingDir, caseName, caseExt,
                                   nBins, logBins, binProp,
                                   regionName=regionName,
                                   dataRegionX=dataRegionX,
                                   dataRegionY=dataRegionY,
                                   recircCenter=recircDefinedRegion,
                                   autoRegion=autoRegion,
                                   includePillar=includePillar,
                                   halfRegion=halfRegion)
print(outFile)
for fileName in fileList:
    if re.match(filePat, fileName):
        print(fileName)
        # Check for fileName already in metaData, skip if so
        if testMode: # Test mode requires you specify exactly which file you wish to test
            data = dh.dataLoader(caseName+caseExt, type=caseExt[2:-1])
        else:
            data = dh.dataLoader(fileName, type=caseExt[2:-1])
        params = dh.extractParams(fileName, nPil, caseExt=caseExt[2:-1])
        if autoRegion:
            dataRegionX, dataRegionY = dh.pillarGapCalculation(params['r1'], params['r2'], params['d'],
                                                            includePillar=includePillar, halfRegion=halfRegion)
        data = dh.subSelectData(data, xRange=dataRegionX, yRange=dataRegionY)
        if testMode & plotData:
            ax1 = plotDataSet(data, "raw")
        if recircDefinedRegion:
            #data = selectRecircZoneNegVel(data, params['r1'], params['r2'], params['d'], gridSize=10)
            data = dh.subSelectData(data,xRange=[250, 500]) # Subselect only half the recirculation area
            data, coords = calcRecircBoundsFluxInt(data, gridSize=1)
            if data.empty:  # SKIP FILES THAT HAVE NO RECIRCULATION
                print("EMPTY FILE")
                continue
        # else:
        # This uses negative velocity criteria to find the edge. Likely not worth using.
        #     data, params['xEdge'] = selectRecircZoneBasic(data, params['r1'], params['r2'], params['d'], includePillar)
        if testMode & plotData:
            ax2 = plotDataSet(data, "Subselected")
        params['recircVol'] = data.eleVol.sum()
        params['dP'], params['q'], params['l'] = calcFlowPress(data, params)
        if estimateRecircCenter:
            params['recircCenter'] = centerPointEstimation(data, r1=params['r1'], r2=params['r2'], d=params['d'], useMid=useMid)
        if testMode & plotData:
            plotPoints(ax2, [params['recircCenter'][0], params['recircCenter'][1], 110])
        if estimateRecircCenter:
            params['totalRecircFlux'], params['posFlux'], params['negFlux'] = \
                estimateRecircFlux(data, params['recircCenter'], params['recircVol'])
            params['posMRT'] = params['recircVol']/params['posFlux']
            params['negMRT'] = params['recircVol']/abs(params['negFlux'])
        else:
            params['totalRecircFlux'] = 0
            params['posFlux'] = 0
            params['negFlux'] = 0
            params['posMRT'] = 0
            params['negMRT'] = 0
        params['estRT'] = params['d']*10**-6/data.velMag.mean()
        params['fileName'] = fileName
        if vortAng:
            data.loc[:, 'angle'] = calcVortVelAngle(data, 'u', 'v', 'w',
                                                    'vortX', 'vortY', 'vortZ')
        elementVol = data.eleVol.values
        params['totalVol'] = np.sum(elementVol)
        params['velChar'] = params['Re']*nu/500E-6  #Assume 500 um channel width
        params['nu'] = nu
        if calcChem:
            kVal = data.loc[data.index[0], 'k']
            dCdt = data.h2o2.values*data.tcpo.values*kVal
            data.loc[:, 'dCdt'] = dCdt
            params['totalProd'] = np.sum(np.multiply(data.cProduct.values, elementVol))
            params['totalTCPO'] = np.sum(np.multiply(data.tcpo.values, elementVol))
            params['totalH2O2'] = np.sum(np.multiply(data.h2o2.values, elementVol))
            params['dCdtMax'] = np.max(dCdt)
            params['dCdtSum'] = np.sum(dCdt)
            params['dilutionTCPO'], params['reactorTCPO'] = calcDilutionIndex(data, 'tcpo')
            data.loc[:, 'constC'] = data.tcpo.values+data.cProduct.values # Conservative component
            params['dilutionConserv'], params['reactorConserv'] = calcDilutionIndex(data, 'constC')
            dCdtNorm = data.h2o2.values*data.tcpo.values/((params['c']/2)**2)  # Max rate is defined by well mixed, which is going to be cA*cB/(c0/2)**2
            data.loc[:, 'dCdtNorm'] = dCdtNorm
            params['maxRef'] = maxValue
            dCdtMaxNorm = data.h2o2.values*data.tcpo.values/maxValue
            data.loc[:, 'dCdtMaxNorm'] = dCdtMaxNorm
            params['conservative'] = np.sum(data.constC.values*data.eleVol.values)
            params['DH2O2'] = diff
            params['Pe'] = params['velChar']*params['r1']*2E-6/diff
            params['DaDiff'] = params['k']*params['c']/1000*(2E-6*params['r1'])**2/diff
            params['DaAdv'] = params['k']*params['c']/1000*2E-6*params['r1']/params['velChar']
        if binProp:
            # Calculate stats of mean and std dev based solely on prop weighted by vol
            params['volWeightedMean'] = np.average(data.loc[:, binProp].values, weights=elementVol)
            var = data.loc[:,binProp].values-params['volWeightedMean']
            params['volWeightedStd']  = np.sqrt(np.average(var**2, weights=elementVol))
            normFreq, valMean, valBin = \
                producePDF(data, nBins=nBins, logBin=logBins, prop=binProp)
            pdfData = {'normFreq': normFreq, 'valMean': valMean,
                       'leftBin': valBin[:-1], 'rightBin': valBin[1:]}
            # Write PDF mean and std to params file, will need to import PDFstats function from elsewhere
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
        if testMode:
            if recircDefinedRegion:
                f4, ax4 = plt.subplots(subplot_kw={"projection":"3d"})
                # ax4.scatter(boundCoords[:,0], boundCoords[:,1], boundCoords[:,2]
                boundCoords = coords[coords[:, 0] != 0]
                x = boundCoords[:, 0]
                y = boundCoords[:, 1]
                z = boundCoords[:, 2]
                ax4.scatter(x, y, z)
                ax4.set_xlabel('x coord')
                ax4.set_ylabel('y coord')
                ax4.set_zlabel('z coord')
                ax4.set_title('Boundary Points')
            plt.ion()
            plt.show()
            print('BREAK')
            break

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
# plt.legend(loc=0)
# plt.ylabel('Pressure Difference (Pa)')
# plt.xlabel('Flow rate (m^3/s)')
# plt.title('Flow rate vs Pressure')
# plt.savefig(caseName+'_dPvsQ.png', dpi=600)
