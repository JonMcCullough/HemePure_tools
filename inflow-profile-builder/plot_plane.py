import numpy as np

@np.vectorize
def calcz(x,y,a,b,c,d):
  return (-a*x - b*y - d) / c

def plot_xyz_plane_3views(tdf):
  """ Plot 3 views of a plane and save it to pngs.
  """
  plot = tdf.plot(x="x",y="y")
  fig = plot.get_figure()
  fig.savefig("xy.png")

  plot = tdf.plot(x="x",y="z")
  fig = plot.get_figure()
  fig.savefig("xz.png")


  plot = tdf.plot(x="y",y="z")
  fig = plot.get_figure()
  fig.savefig("yz.png")

def plot_fitting_plane(ax,tdf,a,b,c,d):
  """ plot a fitting plane.
  """
  xx = np.linspace(tdf["x"].min(), tdf["x"].max(), 1)
  yy = np.linspace(tdf["y"].min(), tdf["y"].max(), 1)
  zz = calcz(xx[:,None], yy[None,:],a,b,c,d)

  xx, yy = np.meshgrid(range(tdf["x"].min(), tdf["x"].max()), range(tdf["y"].min(), tdf["y"].max(),))
  # calculate corresponding z
  zz = (-a * xx - b * yy - d) * 1. / c

  ax.plot_surface(xx, yy, zz, alpha=0.8, color=[0,1,0])

