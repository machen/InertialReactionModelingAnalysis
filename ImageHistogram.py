from PIL import Image
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import re
import seaborn as sns
import DataHelper as dh


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


def genPDF(image, bins):
    pdf, edges = np.histogram(image, density=True, bins=bins)
    edgeVal = (edges[:-1]+edges[1:])/2
    return pdf, edgeVal, edges[:-1], edges[1:]


def pdfStats(data):
    dx = data.rightBin-data.leftBin
    eVal = np.sum(data.valMean*data.normFreq*dx)
    varVal = np.sum(data.normFreq*(data.valMean-eVal)**2*dx)
    return eVal, varVal


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
    params = dh.extractExptParams(file)
    img = Image.open(file)
    imageDictInv = {v: k for k, v in imageDict.items()}
    channelPat = re.compile('C=(\d)')
    channel = int(re.search(channelPat, file).group(1))
    channelName = imageDictInv[channel]
    # Sums on large images can cause overflow issues leading to negative values
    # Forcing type to uInt64 ensures pandas correctly reads the value when converting to dataframe for export.
    # TODO: It may be smarter to normalize this to the image bitdepth so that we deal with floats instead.
    data = np.array(img, dtype='uint64')
    if maxNorm:  # If using maximum normalization, all values scaled to largest value
        if maxVal:
            data = data/maxVal
        else:
            data = data/np.max(data)
    data = dh.subSelectExptData(data, xRange=xRange, yRange=yRange)
    params['fileName'] = os.path.splitext(file)[0]
    params['maxInt'] = np.max(data)
    params['meanInt'] = np.mean(data)
    params['stdInt'] = np.std(data)
    params['sumInt'] = np.sum(data)
    params['channel'] = channelName
    params['caseExt'] = f'.{channelName}_hist'
    dataPdf, dataVal, dataLeft, dataRight = genPDF(data, bins)
    dataDict = {'normFreq': dataPdf, 'valMean': dataVal,
                'leftBin': dataLeft, 'rightBin': dataRight}
    dataDF = pd.DataFrame(dataDict)
    dataDF.to_csv(outFile+params['fileName']+params['caseExt'])
    dataImage = Image.fromarray(data.astype('uint16'))
    dataImage.save(outFile+file[:-4]+'.tiff')
    return params


workingDir = "..\\..\\Experiments\\2023-5-5-ChemilumRandMP\\MPD3B_C3\\StitchSplit\\"
#workingDir = "C:\\Users\\mache\\Google Drive Workspace\\2022-1-15-Chemilum 100 um\\100 um Gap\\SplitImgs\\"
sequenceFile = "..\\RawData\\Sequence1.txt"
os.chdir(workingDir)
#"MP_P3_D1_3c_100q.nd2 - MP_P3_D1_3c_100q.nd2 (series 1) - C=0"
#filePat = re.compile('.*(series 2).*\\.tif')
filePat = re.compile('.*\\.tif')

# Folder of names to further restrict analysis. Uses a text list of filenames
# Set to none to not use
if sequenceFile:
    with open(sequenceFile, 'r') as filterFile:
        filterList = [item.split()[0] for item in filterFile]


bins = 50
# Remember that its is supposed to be the frame of the image
xRange =  [2004, 2331]
yRange = [1880, 2100]
maxNorm = False
# Set to none to use max observed in image, otherwise use well mixed value
maxVal = 920
regionName = "Batch1 Mid Pore B"

fileList = os.listdir()
# Links identifier to stack position, also calls what images will be binned
imageDict = {'bright': 0, 'dark': 1, 'fluor': 2}
outFile = genOutputFolderAndParams(workingDir, filePat, bins, maxNorm, maxVal,
                                   regionName=regionName,
                                   imageDict=imageDict, dataRegionX=xRange,
                                   dataRegionY=yRange)
metaData = []

for file in os.listdir():
    # If you specify a secondary list to filter the items by, this will kick out anything not in that list
    if filterList:
        if any(item in file for item in filterList):
            if re.match(filePat, file):
                print(file)
                params = produceSinglePDF(file, imageDict, outFile, maxNorm, maxVal=maxVal,
                                        bins=100, xRange=xRange, yRange=yRange)
                metaData.append(params)
    else:
        if re.match(filePat, file):
            print(file)
            params = produceSinglePDF(file, imageDict, outFile, maxNorm, maxVal=maxVal,
                                    bins=100, xRange=xRange, yRange=yRange)
            metaData.append(params)
outDF = pd.DataFrame.from_records(metaData, coerce_float=True)
outDF.to_csv(outFile+"_meta.csv")
