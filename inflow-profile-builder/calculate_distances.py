import numpy as np
import pandas as pd
import sys

def calculate_distances_from_borders(tdf, border_points):
  """
  taken a mesh and a set of border points, calculate the distance to the nearest border point.
  """
  tdfa = tdf.to_numpy()
  tdfb = border_points.to_numpy()

  distances = np.zeros((len(tdfa),1))
  distances.fill(1000000.0)

  for i in range(0, len(tdfa)): #JM was xrange
    for j in range(0, len(tdfb)): #JM was xrange
      # 1-liner to calculate distance between two points.
      dist = np.linalg.norm(tdfa[i]-tdfb[j])
      if dist < distances[i]:
        distances[i] = dist
        if dist == 0.0:
          break

  return distances

def calculate_distances_from_borders_per_inlet(tdf, border_points, numInlets):
  """
  taken a mesh and a set of border points, calculate the distance to the nearest border point.
  """
  distances = np.zeros((len(tdf.index),2))
  distances.fill(1000000.0)
  
  iLast = 0
  iCount = 0

  for iNum in range(0,numInlets):
      print(iNum)
      tdfa = tdf[tdf["iolet_number"]==iNum]
      tdfa = tdfa.to_numpy()
     
      tdfb = border_points[border_points["iolet_number"]==iNum]
      tdfb = tdfb.to_numpy()

      for i in range(0, len(tdfa)): #JM was xrange
        for j in range(0, len(tdfb)): #JM was xrange
          # 1-liner to calculate distance between two points.
          dist = np.linalg.norm(tdfa[i]-tdfb[j])
          if dist < distances[iCount,0]:
            distances[iCount,0] = dist
 
            if dist == 0.0:
              break
        iCount = iCount + 1
      
      if iLast<iCount:
          distances[iLast:iCount,1] = max(distances[iLast:iCount,0])
      else:
          print("Size Error on outlet ",iNum)

      iLast = iCount

  return distances

