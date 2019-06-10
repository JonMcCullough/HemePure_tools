import numpy as np

def calc_dist_to_plane(x, y, z, a, b, c, d):
    """
    Projects the points with coordinates x, y, z onto the plane
    defined by a*x + b*y + c*z + d = 0
    """
    
    distance = (a*x + b*y + c*z + d) / np.sqrt(a*a + b*b + c*c)
    return distance



def collapse_to_plane(tdf,a,b,c,d):
  """ Take a point cloud, and collapse it to the plane specified.
  """

  px = []
  py = []
  pz = []

  for i in xrange(0,len(tdf["x"])):
    dist = calc_dist_to_plane(tdf["x"][i], tdf["y"][i], tdf["z"][i], a, b, c, d)
    norm = np.linalg.norm([a,b,c])
    px.append(tdf["x"][i] - (dist * a/norm))
    py.append(tdf["y"][i] - (dist * b/norm))
    pz.append(tdf["z"][i] - (dist * c/norm))

  tdf["px"] = px
  tdf["py"] = py
  tdf["pz"] = pz

  return tdf
