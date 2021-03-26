from PIL import Image
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import re

""" Should read in images (which should be in the order brightfield, dark, GFP)
and then produce a PDF of that image.


Things I need to do:

-To prep for analysis, I need to also do the statistical calculations I
used in the sim analysis suite
-The images have noticeable noise, I should provide ameans to do do background subtraction
-I should consider setting a threshold value which represents "non channel" space.
Any value less than that should get thrown out for the PDF analysis. Would be nice to show it as well

-If I'm feeling really fancy I could use the tracer to give me the indicies for areas that are 0

"""


def genPDF(image, bins):
    pdf, edges = np.histogram(image, density=True, bins=bins)
    edgeVal = (edges[:-1]+edges[1:])/2
    return pdf, edgeVal, edges[:-1], edges[1:]


def pdfStats(data):
    dx = data.rightBin-data.leftBin
    eVal = np.sum(data.valMean*data.normFreq*dx)
    varVal = np.sum(data.normFreq*(data.valMean-eVal)**2*dx)
    return eVal, varVal


def subSelectData(data, xRange=None, yRange=None):
    """Assumes that each of the inputs to the function is a tuple containing
     max and min values of x, y, and z that we wish to include.

     Assumes that we're working with images. You should 100% pre run without any
     subselection of data to determine the correct area"""
    if xRange:
        data = data[:, xRange[0]:xRange[1]]
    if yRange:
        data = data[yRange[0]:yRange[1], :]
    return data


def genOutputFolderAndParams(dataDir, case, nBins, maxNorm,
                             imageDict=None,
                             regionName='Result', dataRegionX=None,
                             dataRegionY=None, dataRegionZ=None):
    # Creates a new folder and populates it with key binning info
    if maxNorm:
        outputPath = '..\\{} max norm {} bins\\'.format(regionName, nBins)
    else:
        outputPath = '..\\{} {} bins\\'.format(regionName, nBins)
    outputFile = outputPath+'Parameters.txt'
    if not os.path.isdir(outputPath):
        os.mkdir(outputPath)
    with open(outputFile, "w") as outFile:
        outFile.write("Data directory: {}\n".format(dataDir))
        outFile.write("Case extension: {}\n".format(case))
        outFile.write("Number of bins: {}\n".format(nBins))
        outFile.write("Image key: {}\n".format(imageDict))
        outFile.write("Maximum Normalization: {}\n".format(maxNorm))
        outFile.write("X Region: {}\n".format(dataRegionX))
        outFile.write("Y Region: {}\n".format(dataRegionY))
        outFile.write("Z Region: {}\n".format(dataRegionZ))
    return outputPath


def produceSinglePDF(file, imageDict, outFile, maxNorm, maxVal=None, bins=100,
                     xRange=None, yRange=None):
    # Single image. No background subtraction for the moment.
    img = Image.open(file)
    imageDictInv = {v: k for k, v in imageDict.items()}
    channelPat = re.compile('C=(\d)')
    channel = int(re.search(channelPat, file).group(1))
    channelName = imageDictInv[channel]
    data = np.array(img)
    if maxNorm:  # If using maximum normalization, all values scaled to largest value
        if maxVal:
            data = data/maxVal
        else:
            data = data/np.max(data)
    data = subSelectData(data, xRange=xRange, yRange=yRange)
    dataPdf, dataVal, dataLeft, dataRight = genPDF(data, bins)
    dataDict = {'normFreq': dataPdf, 'valMean': dataVal,
                'leftBin': dataLeft, 'rightBin': dataRight}
    dataDF = pd.DataFrame(dataDict)
    dataDF.to_csv(outFile+file[:-4]+'_{}_hist.csv'.format(channelName))
    dataImage = Image.fromarray(data)
    dataImage.save(outFile+file[:-4]+'.tiff')
    return


workingDir = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\Mar22_2021-Chemilum\\2PD-1_P4_A1 - 100 um gap\\BackgroundSub\\"
os.chdir(workingDir)
filePat = re.compile('.*\.tif')
bins = 50
xRange = [800, 1140]  # [800, 1300]  # Should be matrix indices for the given image, you must update this
yRange = [920, 1070]  # [900, 1500]  # Should be matrix indices for the given image
maxNorm = False
maxVal = 1748.3  # Set to none to use max observed in image. Value will depend on specific days mix
regionName = "Pillar gap"
bgFile = 'NoDevice.tif'

fileList = os.listdir()
# Links identifier to stack position, also calls what images will be binned
imageDict = {'bright': 0, 'dark': 1, 'fluor': 2}
outFile = genOutputFolderAndParams(workingDir, filePat, bins, maxNorm,
                                   regionName=regionName,
                                   imageDict=imageDict, dataRegionX=xRange,
                                   dataRegionY=yRange)

for file in os.listdir():
    if re.match(filePat, file):
        print(file)
        produceSinglePDF(file, imageDict, outFile, maxNorm, maxVal=maxVal,
                         bins=100, xRange=xRange, yRange=yRange)
        # img = Image.open(file)
        # img.seek(0)
        # bright = np.array(img)
        # brightPdf, brightVal, brightLeft, brightRight = genPDF(bright, bins)
        # brightDict = {'normFreq': brightPdf, 'meanVal': brightVal,
        #               'leftBin': brightLeft, 'rightBin': brightRight}
        # brightDF = pd.DataFrame(brightDict)
        # brightDF.to_csv(workingDir+file[:-4]+'_brightfield_histogram.csv')
        # img.seek(1)
        # dark = np.array(img)
        # darkPdf, darkVal, darkLeft, darkRight = genPDF(dark, bins)
        # darkDict = {'normFreq': darkPdf, 'meanVal': darkVal,
        #             'leftBin': darkLeft, 'rightBin': darkRight}
        # darkDF = pd.DataFrame(darkDict)
        # darkDF.to_csv(workingDir+file[:-4]+'_dark_histogram.csv')
        # img.seek(2)
        # fluor = np.array(img)
        # fluorPdf, fluorVal, fluorLeft, fluorRight = genPDF(fluor, bins)
        # fluorDict = {'normFreq': fluorPdf, 'meanVal': fluorVal,
        #              'leftBin': fluorLeft, 'rightBin': fluorRight}
        # fluorDF = pd.DataFrame(fluorDict)
        # fluorDF.to_csv(workingDir+file[:-4]+'_fluorescence_histogram.csv')
