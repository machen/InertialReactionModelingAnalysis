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
-Also, probably useful to window my data (may be possible to eliminate background?)
Could avoid with log PDF


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


workingDir = "G:\\My Drive\\Postdoctoral work\\Inertial flow study\\Experiments\\2PillarD-1_P2_A2\\DataProcessing\\"
os.chdir(workingDir)
filePat = re.compile('.*\.tif')
bins = 100

fileList = os.listdir()


for file in os.listdir():
    if re.match(filePat, file):
        print(file)
        img = Image.open(file)
        img.seek(0)
        bright = np.array(img)
        brightPdf, brightVal, brightLeft, brightRight = genPDF(bright, bins)
        brightDict = {'normFreq': brightPdf, 'meanVal': brightVal,
                      'leftBin': brightLeft, 'rightBin': brightRight}
        brightDF = pd.DataFrame(brightDict)
        brightDF.to_csv(workingDir+file[:-4]+'_brightfield_histogram.csv')
        img.seek(1)
        dark = np.array(img)
        darkPdf, darkVal, darkLeft, darkRight = genPDF(dark, bins)
        darkDict = {'normFreq': darkPdf, 'meanVal': darkVal,
                    'leftBin': darkLeft, 'rightBin': darkRight}
        darkDF = pd.DataFrame(darkDict)
        darkDF.to_csv(workingDir+file[:-4]+'_dark_histogram.csv')
        img.seek(2)
        fluor = np.array(img)
        fluorPdf, fluorVal, fluorLeft, fluorRight = genPDF(fluor, bins)
        fluorDict = {'normFreq': fluorPdf, 'meanVal': fluorVal,
                     'leftBin': fluorLeft, 'rightBin': fluorRight}
        fluorDF = pd.DataFrame(fluorDict)
        fluorDF.to_csv(workingDir+file[:-4]+'_fluorescence_histogram.csv')
