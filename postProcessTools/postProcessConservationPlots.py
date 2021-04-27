import matplotlib
matplotlib.use('Agg')

import os,sys

from postProcessPlotSingleFunction import hemeCollateFunction

import numpy as np
import matplotlib.pyplot as plt
import tikzplotlib

def execute(command):
        print("Executing: " + command)
        r = os.system(command)
        if r != 0:
                sys.exit("Command failed.")

#CHANGE for Specific case
stepSize = 10000
DX = 0.000025

#Location (LU) and area (PU) of IOlets
FirstIn = np.array([[1352.64,538.63,10766.2]])
FirstInArea = np.array([3.84949e-6])

FirstOut = np.array([[30.0271,567.329,13.128],[432.466,1476.36,1637.82],[596.032,997.209,4475.71]])
FirstOutArea = np.array([1.46699e-06,2.00232e-06,9.32602e-07])

if len(sys.argv) != 2:
    sys.exit("usage: python3 hemepreprocName Folder")

dataSet = sys.argv[1]



data = np.load(dataSet+'.npz')

FirstInV = data['a']
FirstOutV = data['b']


names=dataSet.split('_')
datalength = len(FirstInV[0][1::3])

######################################################
plt.figure(4)
plt.title("Mean Flow Rate")
for io in range(np.shape(FirstIn)[0]):
    plt.plot(np.linspace(0,datalength/40,datalength),FirstInV[io][1::3]*FirstInArea[io],'-', label='Artery Inlet')

TotalQout = np.zeros(len(FirstOutV[0][1::3]))
for io in range(np.shape(FirstOut)[0]):
    TotalQout += FirstOutV[io][1::3]*FirstOutArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),TotalQout,'-', label='Total Artery Outlet')
plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_InOutArtery_Total_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_InOutArtery_Total.png")
tikzplotlib.save("fig_"+names[0]+"_InOutArtery_Total_"+names[1]+".tex")

######################################################
plt.figure(7)
plt.title("Mean Flow Rate")

TotalQout = np.zeros(len(FirstOutV[0][1::3]))
for io in range(np.shape(FirstOut)[0]):
    TotalQout += FirstOutV[io][1::3]*FirstOutArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),np.cumsum(TotalQout),'-', label='Total Artery Outlet')

plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_Beds_Cum_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_Beds_Cum.png")
tikzplotlib.save("fig_"+names[0]+"_Beds_Cum_"+names[1]+".tex")

