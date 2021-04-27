import matplotlib
matplotlib.use('Agg')

from postProcessPlotSingleFunction import hemeCollateFunction

import numpy as np
import matplotlib.pyplot as plt

#CHANGE THESE FOR SPECIFIC DATA
stepSize = 100
DX = 0.0002

#Location (LU) and area (PU) of IOlets 
FirstIn = np.array([[23,8,3]])
FirstInArea = np.array([3.14e-6])

FirstOut = np.array([[5.5,8,43.4],[23,8,53.5],[40.5,8,43.4]])
FirstOutArea = np.array([3.14e-6,3.14e-6,3.14e-6])

FirstInP, FirstInV = hemeCollateFunction("results_0/inletFirst.txt", stepSize, DX, FirstIn)
FirstOutP, FirstOutV = hemeCollateFunction("results_0/outletFirst.txt", stepSize, DX, FirstOut)

SimSteps = np.linspace(0,stepSize*len(FirstInV[0][1::3]),len(FirstInV[0][1::3]))

######################################################
plt.figure(1)
plt.title("Flow rate")
for io in range(np.shape(FirstIn)[0]):
    plt.plot(SimSteps,FirstInV[io][1::3]*FirstInArea[io],'-', label='Total Inlet')

plt.xlabel("Simulation Steps")
plt.ylabel("Volumetric Flow Rate [m3/s]")
plt.legend()
plt.savefig("fig_InOutFlow.png")

######################################################
plt.figure(2)
plt.suptitle("Artery Outlet")
plt.subplot(2,1,1)
for io in range(np.shape(FirstOut)[0]):
    plt.plot(FirstOutV[io][1::3],'-', label='IOlet #'+str(io))
#plt.title("Velocity")
#plt.xlabel("Simulation Steps")
plt.ylabel("Velocity [Units]")
plt.legend()

plt.subplot(2,1,2)
for io in range(np.shape(FirstOut)[0]):
    plt.plot(FirstOutP[io][1::3],'-', label='IOlet #'+str(io))
#plt.title("Pressure")
plt.xlabel("Simulation Steps")
plt.ylabel("Pressure [Units]")
plt.legend()

plt.savefig("fig_FirstOutlet.png")

#plt.show()
