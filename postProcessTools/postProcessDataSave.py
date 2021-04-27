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

if len(sys.argv) != 2:
    sys.exit("usage: python3 hemepreprocName Folder")

dataSet = sys.argv[1]

FirstInP, FirstInV = hemeCollateFunction(dataSet+"/inletFirst.txt", stepSize, DX, FirstIn)
FirstOutP, FirstOutV = hemeCollateFunction(dataSet+"/outletFirst.txt", stepSize, DX, FirstOut)

DataSite = dataSet+".npz"

np.savez(DataSite,a=FirstInV,b=FirstOutV,c=SecondInV,d=SecondOutV)

