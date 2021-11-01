import numpy as np
import pandas as pd

def count_borders(border_ids, links, i):
  borders = 0
  for i in links[i]:
    if i in border_ids:
      borders += 1
  return borders



def detect_border_points_in_plane(tdf, links, threshold=10):

  # list of IDs of border points found.
  border_ids = []

  xs = pd.DataFrame(tdf, columns=['px', 'py', 'pz']).to_numpy()

  # initial border selection
  border_points = np.zeros((0,3))
  for i in range(0, len(xs)): #JM was xrange
    if len(links[i]) < threshold:
      border_ids.append(i)

  # neighbour check, we prune points that are not close to other border points
  #tmp_border_ids = []
  #for i in xrange(0, len(border_ids)):
  #  if count_borders(border_ids, links, i) > 0:
  #    tmp_border_ids.append(i)

  #border_ids = tmp_border_ids

  # lastly, we populate the border_points array with the final list of border ids.
  for i in range(0, len(border_ids)): #JM was xrange
    border_points = np.vstack((border_points, xs[border_ids[i]]))

  return border_points, border_ids
