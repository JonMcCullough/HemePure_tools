import numpy as np
import sys

path = sys.argv[1]

print("starting to process", path)
f = open(path, "r")

wtSum = 0
wtCount = 0
wtMax = 0
wtMin = 10
for i in f:  
    i = np.float_(i.split())

    wtSum = wtSum + i[3]
    wtCount = wtCount + 1

    if i[3] < wtMin:
        wtMin = i[3]

    if i[3] > wtMax:
        wtMax = i[3]

print("For ", wtCount, " points:")
print("Max = ", wtMax, ", Min = ", wtMin, ", Ave = ", wtSum/wtCount)

