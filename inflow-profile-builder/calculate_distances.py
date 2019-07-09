import numpy as np
import pandas as pd
import sys

def calculate_distances_from_borders(tdf, border_points):
  """
  taken a mesh and a set of border points, calculate the distance to the nearest border point.
  """
  tdfa = tdf.as_matrix()
  tdfb = border_points.as_matrix()

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

