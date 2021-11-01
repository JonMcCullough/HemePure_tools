import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import scipy.optimize
import functools
from scipy.optimize import leastsq


def numpyify(points):
  return points.to_numpy()

def f_min(X,p):
    plane_xyz = p[0:3]
    distance = (plane_xyz*X.T).sum(axis=1) + p[3]
    return distance / np.linalg.norm(plane_xyz)

def residuals(params, signal, X):
    return f_min(X, params)

def fit_plane(amb_points, plot=False):

  # returns a,b,c,d in ax+by+cz+d=0. a,b,c are also the normal.
  points = numpyify(amb_points)

  pointsT = points.transpose()
  # Inital guess of the plane
  diff = points[0] - points[-1]

  p0 = np.array(([diff[0], diff[1], diff[2], 1.]))

  sol = leastsq(residuals, p0, args=(None, pointsT))[0]

  #print "Solution: ", sol
  #print "Old Error: ", (f_min(pointsT, p0)**2).sum()
  #print "New Error: ", (f_min(pointsT, sol)**2).sum()

  return sol

