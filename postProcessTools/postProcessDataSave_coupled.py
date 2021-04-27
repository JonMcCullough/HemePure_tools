import matplotlib
matplotlib.use('Agg')

import os,sys

from postProcessPlotSingleFunction import hemeCollateFunction

import numpy as np
import matplotlib.pyplot as plt

def execute(command):
        print("Executing: " + command)
        r = os.system(command)
        if r != 0:
                sys.exit("Command failed.")

#CHANGE for specific case
stepSize = 10000
DX = 0.000025

#Location of IOlets (LU)
FirstIn = np.array([[1352.64,538.63,10766.2]])

FirstOut = np.array([[30.0271,567.329,13.128],[432.466,1476.36,1637.82],[596.032,997.209,4475.71]])

SecondIn = np.array([[48.5779,1246.65,65.275],[72.1172,1466.24,289.857],[99.6054,1906.81,203.998],[136.753,796.972,126.79],[231.314,2143.19,105.563],[525.19,647.567,67.4014],[606.887,1979.7,48.1466],[615.394,2801.43,86.8552],[661.903,1290.98,7.4712],[736.591,2316.31,33.8911],[739.486,1901.7,70.8137],[903.65,2464.36,5559.7],[1044.49,2372.05,229.391],[1068.12,2880.39,776.374],[1344.01,3406.65,1027.61],[1516.86,3575.37,814.952],[1636.03,3499.42,1255.14]])

SecondOut = np.array([[253.017,1243.39,11673],[289.949,1486.25,11676.4],[295.58,851.12,11677.3],[765.8,258.572,11677.8],[1067.79,2289.71,11658.8],[1177.44,84.8546,11676.9],[2133.53,1416.29,11677.4],[2151.22,1858.7,11666.2],[2184.71,1015.9,11677.6],[2194.39,15.9046,11676.6],[2230.55,1512.43,11668.5],[2649.39,625.144,11676.9],[2706.98,1443.99,11674.9],[2735.5,971.116,11663.6],[2996.37,505.185,11674.8]])


if len(sys.argv) != 2:
    sys.exit("usage: python3 hemepreprocName Folder")

dataSet = sys.argv[1]

FirstInP, FirstInV = hemeCollateFunction(dataSet+"/inletFirst.txt", stepSize, DX, FirstIn)
FirstOutP, FirstOutV = hemeCollateFunction(dataSet+"/outletFirst.txt", stepSize, DX, FirstOut)

SecondOutP, SecondOutV = hemeCollateFunction(dataSet+"/outletSecond.txt", stepSize, DX, SecondOut)
SecondInP, SecondInV = hemeCollateFunction(dataSet+"/inletSecond.txt", stepSize, DX, SecondIn)

DataSite = dataSet+".npz"

np.savez(DataSite,a=FirstInV,b=FirstOutV,c=SecondInV,d=SecondOutV)

