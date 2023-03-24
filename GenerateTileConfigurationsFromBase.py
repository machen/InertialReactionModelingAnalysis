
""" Script is used to generate tile configurations that can
be used to stitch together images. This relies on stitching first
one image (say a brightfield image) and getting a TileConfiguration
file from it. This then substitutes the appropriate parts of the
tile configuration and writes out a new one which is named after the
file it should be applied to.


"""
import os
from string import Template

baseDir = "..\\..\\Experiments\\2023-3-21-Chemilum\\MPD1_D4_Batch2\\"
imageDir = "RawData\\" # Directory of images that gives the base image names. Should be the raw data folder
writeLoc = "SplitImg\\" # Where do you want the template files to be put?
ext = ".nd2" # Extension of the raw data
os.chdir(baseDir)

# This should be a registered tile configuration produced by stitching a brightfield image.
tempLoc = "SplitImg\\TileConfigurationTemplate.txt"
with open(tempLoc, 'r') as tempFile:
    tempStr = tempFile.read()
    template = Template(tempStr)

for imageName in os.listdir(imageDir):
    if imageName.endswith(ext):
        for i in range(0,3):
            outFileName = writeLoc+imageName+"C={}".format(i)+"TileConfig.txt"
            with open(outFileName, 'x') as outFile:
                outFile.write(template.safe_substitute(baseFile=imageName,
                              channelNo=i))
    else:
        continue