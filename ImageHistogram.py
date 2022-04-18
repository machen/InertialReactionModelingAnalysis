from PIL import Image
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import re
import seaborn as sns
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
plt.rcParams['svg.fonttype'] = 'none'
sns.set_context('talk')


def extractParams(fileName):
    # Produces a dictionary of experimental parameters
    cPat = re.compile('(\d+\.?\d*)c')
    qPat = re.compile('(\d*\.?\d*)q')
    repPat = re.compile('(\d{0,3}).nd2')
    qVal = re.search(qPat, fileName).group(1)
    repVal = re.search(repPat, fileName).group(1)
    cVal = re.search(cPat, fileName).group(1)
    if not repVal:
        repVal = 0
    res = {'q': float(qVal), 'replicate': int(repVal), 'c': float(cVal)}
    return res


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


def genOutputFolderAndParams(dataDir, case, nBins, maxNorm, maxVal,
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
        outFile.write("Maximum Normalization: {}, Value: {}\n".format(maxNorm, maxVal))
        outFile.write("X Region: {}\n".format(dataRegionX))
        outFile.write("Y Region: {}\n".format(dataRegionY))
        outFile.write("Z Region: {}\n".format(dataRegionZ))
    return outputPath


def produceSinglePDF(file, imageDict, outFile, maxNorm, maxVal=None, bins=100,
                     xRange=None, yRange=None):
    # Single image. No background subtraction for the moment.
    params = extractParams(file)
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
    params['fileName'] = os.path.splitext(file)[0]
    params['maxInt'] = np.max(data)
    params['meanInt'] = np.mean(data)
    params['stdInt'] = np.std(data)
    params['channel'] = channelName
    params['caseExt'] = f'.{channelName}_hist'
    dataPdf, dataVal, dataLeft, dataRight = genPDF(data, bins)
    dataDict = {'normFreq': dataPdf, 'valMean': dataVal,
                'leftBin': dataLeft, 'rightBin': dataRight}
    dataDF = pd.DataFrame(dataDict)
    dataDF.to_csv(outFile+params['fileName']+params['caseExt'])
    dataImage = Image.fromarray(data)
    dataImage.save(outFile+file[:-4]+'.tiff')
    return params


workingDir = "..\\..\\Experiments\\2022-3-22-MPD2\\MPD2_P1_A3\\StitchImgs\\"
#workingDir = "C:\\Users\\mache\\Google Drive Workspace\\2022-1-15-Chemilum 100 um\\100 um Gap\\SplitImgs\\"
os.chdir(workingDir)
#"MP_P3_D1_3c_100q.nd2 - MP_P3_D1_3c_100q.nd2 (series 1) - C=0"
#filePat = re.compile('.*(series 5).*\.tif')
filePat = re.compile('.*\.tif')
bins = 50
# Remember that tis is supposed to be the frame of the image
xRange = [696, 8760]
yRange = [1032, 1524]
maxNorm = False
# Set to none to use max observed in image, otherwise use well mixed value
maxVal = 1083
regionName = "Raw Masked Pore Throats"

fileList = os.listdir()
# Links identifier to stack position, also calls what images will be binned
imageDict = {'bright': 0, 'dark': 1, 'fluor': 2}
outFile = genOutputFolderAndParams(workingDir, filePat, bins, maxNorm, maxVal,
                                   regionName=regionName,
                                   imageDict=imageDict, dataRegionX=xRange,
                                   dataRegionY=yRange)
metaData = pd.DataFrame([], columns=['fileName', 'c', 'q'])

for file in os.listdir():
    if re.match(filePat, file):
        print(file)
        params = produceSinglePDF(file, imageDict, outFile, maxNorm, maxVal=maxVal,
                                  bins=100, xRange=xRange, yRange=yRange)
        metaData = metaData.append(params, ignore_index=True)
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

metaData.to_csv(outFile+"_meta.csv")
