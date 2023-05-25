from PIL import Image
import numpy as np
import os
import pandas as pd
import DataHelper as dh


workingDir = "..\\..\\Experiments\\2023-4-24-Chemilum-RandPM\\MPD3B_C1\\StitchSplit\\"
os.chdir(workingDir)
# List specifying what images to compare. All flow rates must be unique.

sequenceFile = "..\\RawData\\Batch2Seq1.txt"
with open(sequenceFile, 'r') as filterFile:
    filterList = [item.split()[0] for item in filterFile]

# Use to examine a smaller version of the image
xRange =  [2004, 2331]
yRange = [1880, 2100]
diffRange = [-500, 500]
bitDepth = 16  # Must be set to the image bit depth


regionName = "Batch1 Whole Image"


params = {}
images = {}


# Create a dict with the image parameters and the images. File names are the key values.

#TODO: Needs to consider the fact that FIJI processing adds stuff. Need to figure this better.

for item in filterList:
    params[item] = dh.extractExptParams(item)
    with Image.open(item) as image:
        # Convert image data to float to reduce risk of integer shenanigans
        data = np.array(image)/(2**bitDepth-1)
    images[item] = data

metaData = pd.DataFrame(params)
metaData.sort_values(by='q',inplace=True)
# Need to decide how to do the subtractions, do we just subtract the "nearest" flow rates? I can do this if I have a table sorted by flow rate.