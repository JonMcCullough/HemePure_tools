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

SecondIn = np.array([[48.5779,1246.65,65.275],[72.1172,1466.24,289.857],[99.6054,1906.81,203.998],[136.753,796.972,126.79],[231.314,2143.19,105.563],[525.19,647.567,67.4014],[606.887,1979.7,48.1466],[615.394,2801.43,86.8552],[661.903,1290.98,7.4712],[736.591,2316.31,33.8911],[739.486,1901.7,70.8137],[903.65,2464.36,5559.7],[1044.49,2372.05,229.391],[1068.12,2880.39,776.374],[1344.01,3406.65,1027.61],[1516.86,3575.37,814.952],[1636.03,3499.42,1255.14]])
#SecondInArea = np.array([3.54428e-07,1.25687e-07,1.65423e-07,4.65353e-07,4.13883e-07,6.54601e-07,2.18559e-07,7.61498e-07,1.17853e-07,1.68204e-07,5.6365e-07,8.95763e-07,1.87385e-07,1.33792e-07,7.15095e-08,1.61083e-07,7.39303e-08])#NARROW
SecondInArea = np.array([3.54428e-07,1.25687e-07,1.65423e-07,4.65353e-07,4.13883e-07,6.54601e-07,2.18559e-07,7.61498e-07,1.17853e-07,1.68204e-07,5.6365e-07,1.20339e-05,1.87385e-07,1.33792e-07,7.15095e-08,1.61083e-07,7.39303e-08])#DILATE

SecondOut = np.array([[253.017,1243.39,11673],[289.949,1486.25,11676.4],[295.58,851.12,11677.3],[765.8,258.572,11677.8],[1067.79,2289.71,11658.8],[1177.44,84.8546,11676.9],[2133.53,1416.29,11677.4],[2151.22,1858.7,11666.2],[2184.71,1015.9,11677.6],[2194.39,15.9046,11676.6],[2230.55,1512.43,11668.5],[2649.39,625.144,11676.9],[2706.98,1443.99,11674.9],[2735.5,971.116,11663.6],[2996.37,505.185,11674.8]])
#Original ones - SecondOut = np.array([[258.167,1243.68,11716.4],[297.027,1488.28,11716.4],[302.181,850.661,11716.4],[771.896,259.869,11716.4],[1071.6,2289.37,11716.4],[1185.01,84.7393,11716.4],[2141.7,1417.05,11716.4],[2158.67,1859.51,11716.4],[2191,1017.44,11716.4],[2201.37,16.9951,11716.4],[2238.35,1512.69,11716.4],[2655.69,627.409,11716.4],[2713.34,1442.97,11716.4],[2743.54,971.199,11716.4],[3003.47,506.651,11716.4]])
# Original ones - SecondOutArea = np.array([2.94107e-07,3.40615e-07,9.0167e-07,6.45644e-07,3.07487e-07,5.99456e-07,2.44068e-06,5.8673e-06,1.04336e-06,6.39183e-07,3.51093e-06,7.83413e-06,1.23003e-06,2.53927e-06,4.21734e-07])

#SecondOutArea = np.array([2.74755e-07,3.92763e-07,9.04356e-07,6.10632e-07,4.74039e-07,5.71123e-07,2.44276e-06,7.02101e-06,1.04848e-06,6.15006e-07,4.27224e-06,7.77482e-06,1.27259e-06,3.33818e-06,4.0436e-07])#NARROW
SecondOutArea = np.array([2.74755e-07,3.92763e-07,9.04356e-07,6.10632e-07,4.74039e-07,5.71123e-07,2.44276e-06,3.34791e-05,1.04848e-06,6.15006e-07,4.27224e-06,7.77482e-06,1.27259e-06,3.33818e-06,4.0436e-07])#DILATE

if len(sys.argv) != 2:
    sys.exit("usage: python3 hemepreprocName Folder")

dataSet = sys.argv[1]

BedMap = np.array([[[0],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]]], dtype='object')

data = np.load(dataSet+'.npz')

FirstInV = data['a']
SecondOutV = data['d']
FirstOutV = data['b']
SecondInV = data['c']

names=dataSet.split('_')
datalength = len(FirstInV[0][1::3])

######################################################
#plt.figure(1)
#plt.title("Flow rate")
#for io in range(np.shape(FirstIn)[0]):
#    plt.plot(np.linspace(0,datalength/40,datalength),FirstInV[io][1::3]*FirstInArea[io],'-', label='Inlet #'+str(io))
#for io in range(np.shape(SecondOut)[0]):
#    plt.plot(np.linspace(0,datalength/40,datalength),SecondOutV[io][1::3]*SecondOutArea[io],'-', label='Outlet #'+str(io))
#plt.xlabel("Simulation Steps")
#plt.ylabel("Volumetric Flow Rate [m3/s]")
#plt.legend()
#plt.savefig("fig_InOutFlow.png")
#tikzplotlib.save("fig_InOutFlow.tex")

######################################################
plt.figure(2)
plt.title("Mean Flow Rate")
for io in range(np.shape(FirstIn)[0]):
    plt.plot(np.linspace(0,datalength/40,datalength),FirstInV[io][1::3]*FirstInArea[io],'-', label='Artery Inlet')

TotalQout = np.zeros(len(SecondOutV[0][1::3]))
for io in range(np.shape(SecondOut)[0]):
    TotalQout += SecondOutV[io][1::3]*SecondOutArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),TotalQout,'-', label='Total Vein Outlet')
plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_InOutFlow_Total_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_InOutFlow_Total.png")
tikzplotlib.save("fig_"+names[0]+"_InOutFlow_Total_"+names[1]+".tex")

######################################################
plt.figure(3)
plt.title("Cumulative total Flow Rate")
for io in range(np.shape(FirstIn)[0]):
    plt.plot(np.linspace(0,datalength/40,datalength),np.cumsum(FirstInV[io][1::3]*FirstInArea[io]),'-', label='Artery Inlet')

TotalQout = np.zeros(len(SecondOutV[0][1::3]))
for io in range(np.shape(SecondOut)[0]):
    TotalQout += SecondOutV[io][1::3]*SecondOutArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),np.cumsum(TotalQout),'-', label='Total Vein Outlet')
plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_InOutFlow_Cumulative_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_InOutFlow_Cumulative.png")
tikzplotlib.save("fig_"+names[0]+"_InOutFlow_Cumulative_"+names[1]+".tex")

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
plt.figure(5)
plt.title("Mean Flow Rate")
TotalQin = np.zeros(len(SecondInV[0][1::3]))
for io in range(np.shape(SecondIn)[0]):
    TotalQin += SecondInV[io][1::3]*SecondInArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),TotalQin,'-', label='Total Vein Inlet')

TotalQout = np.zeros(len(SecondOutV[0][1::3]))
for io in range(np.shape(SecondOut)[0]):
    TotalQout += SecondOutV[io][1::3]*SecondOutArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),TotalQout,'-', label='Total Vein Outlet')
plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_InOutVein_Total_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_InOutVein_Total.png")
tikzplotlib.save("fig_"+names[0]+"_InOutVein_Total_"+names[1]+".tex")

######################################################
plt.figure(6)
plt.title("Mean Flow Rate")

TotalQout = np.zeros(len(FirstOutV[0][1::3]))
for io in range(np.shape(FirstOut)[0]):
    TotalQout += FirstOutV[io][1::3]*FirstOutArea[io]
#plt.plot(np.linspace(0,datalength/40,datalength),TotalQout,'-', label='Total Artery Outlet')

TotalQin = np.zeros(len(SecondInV[0][1::3]))
for io in range(np.shape(SecondIn)[0]):
    TotalQin += SecondInV[io][1::3]*SecondInArea[io]
#plt.plot(np.linspace(0,datalength/40,datalength),TotalQin,'-', label='Total Vein Inlet')

plt.plot(np.linspace(0,datalength/40,datalength),TotalQin/TotalQout,'-', label='Ratio In/Out')

plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_Beds_Total_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_Beds_Ratio.png")
tikzplotlib.save("fig_"+names[0]+"_Beds_Total_"+names[1]+".tex")

######################################################
plt.figure(7)
plt.title("Mean Flow Rate")

TotalQout = np.zeros(len(FirstOutV[0][1::3]))
for io in range(np.shape(FirstOut)[0]):
    TotalQout += FirstOutV[io][1::3]*FirstOutArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),np.cumsum(TotalQout),'-', label='Total Artery Outlet')

TotalQin = np.zeros(len(SecondInV[0][1::3]))
for io in range(np.shape(SecondIn)[0]):
    TotalQin += SecondInV[io][1::3]*SecondInArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),np.cumsum(TotalQin),'-', label='Total Vein Inlet')

plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_Beds_Cum_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_Beds_Cum.png")
tikzplotlib.save("fig_"+names[0]+"_Beds_Cum_"+names[1]+".tex")

######################################################
plt.figure(8)
plt.title("Mean Flow Rate")

TotalQout = np.zeros(len(FirstOutV[0][1::3]))
for io in range(np.shape(FirstOut)[0]):
    TotalQout += FirstOutV[io][1::3]*FirstOutArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),TotalQout,'-', label='Total Artery Outlet')

TotalQin = np.zeros(len(SecondInV[0][1::3]))
for io in range(np.shape(SecondIn)[0]):
    TotalQin += SecondInV[io][1::3]*SecondInArea[io]
plt.plot(np.linspace(0,datalength/40,datalength),TotalQin,'-', label='Total Vein Inlet')

#plt.plot(np.linspace(0,datalength/40,datalength),TotalQin/TotalQout,'-', label='Ratio In/Out')

plt.xlabel("Heartbeat Cycle")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
#plt.savefig("fig_"+names[0]+"_Beds_Total_"+names[1]+".png")
plt.savefig("fig_"+dataSet+"_Beds_Total.png")
tikzplotlib.save("fig_"+names[0]+"_Beds_Total_"+names[1]+".tex")

#
#execute("rm -rf " + dataSet + "_figs")
#execute("mkdir " + dataSet + "_figs")
#execute("mv fig_* " + dataSet + "_figs")
#execute("cp fig_* " + dataSet + "_figs")
