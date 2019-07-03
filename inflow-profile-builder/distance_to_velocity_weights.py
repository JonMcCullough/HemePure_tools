import numpy as np

def distance_to_velocity_weights(dist):
  """
  Converts the relative distance to the nearest wall site to a weight vector for the velocity inlets.
  """
  max_dist = dist.max()
  vw = np.zeros((len(dist)))
  for i in range(0, len(dist)): #JM was xrange
    rel_dist = dist[i] / max_dist
    # quadratic function, ax^2+bx+c=0, a=-1, b=2, c=0
    vw[i] = -1.0 * (rel_dist*rel_dist) + 2.0 * rel_dist 
  return vw
