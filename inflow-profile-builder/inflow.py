import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import fit_plane
import plot_plane
import transform_points
import detect_borders
import sys


if __name__ == "__main__":

  #test_data = pd.read_csv("test.csv", delim_whitespace=True)
  test_data = pd.read_csv(sys.argv[1], delim_whitespace=True)
  test_data.columns = ["x","y","z","site_type","iolet_number"]
  adf = pd.DataFrame(test_data, columns=["x","y","z"])


  #plot_plane.plot_xyz_plane_3views(tdf)
  border_points = test_data[test_data['site_type'] == 4] 
  bdf = pd.DataFrame(border_points, columns=["x","y","z"])

  inlet_points = test_data[test_data['site_type'] == 2] 
  idf = pd.DataFrame(inlet_points, columns=["x","y","z"])

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(bdf["x"], bdf["y"], bdf["z"])

  #a,b,c,d = fit_plane.fit_plane(idf)
  #plot_plane.plot_fitting_plane(ax,idf,a,b,c,d)
  #plt.show()
  fig.savefig("plane.png")

  #tdf_collapsed = transform_points.collapse_to_plane(tdf,a,b,c,d)
  #fig = plt.figure()
  #ax = fig.add_subplot(111, projection='3d')
  #ax.scatter(tdf["x"], tdf["y"], tdf["z"])

  plt.title('This should be an exact outline of the inlet plane')

  #import create_mesh_links
  #links = create_mesh_links.create_mesh_links(tdf, dist_limit=3.0)

  #border_points, border_point_ids = detect_borders.detect_border_points_in_plane(tdf_collapsed, links, threshold=14) #13 looked pretty good...
  #ax.scatter(border_points[:,0], border_points[:,1], border_points[:,2])
  ax.scatter(bdf["x"], bdf["y"], bdf["z"])

  #print tdf.loc[[0], ["x","y","z"]]
  #plt.show() #JM Plotting figures won't work on HPC resources

  import calculate_distances
  distances = calculate_distances.calculate_distances_from_borders(adf, bdf)

  import distance_to_velocity_weights
  weights = distance_to_velocity_weights.distance_to_velocity_weights(distances)

  #JM Plotting figures won't work on HPC resources  
  #fig = plt.figure()
  #ax = fig.add_subplot(111, projection='3d')
  #ax.scatter(adf["x"], adf["y"], weights)

  #plt.show()

 
 
  for let in range(0,test_data['iolet_number'].max()+1):
      text_file = open("out" + str(let) + ".weights.txt", "w")
      for i in range(0, len(adf["x"])): #JM was xrange
          if test_data['iolet_number'][i] == let:
              text_file.write("%i %i %i %f\n" % (adf["x"][i], adf["y"][i], adf["z"][i], weights[i]))

      text_file.close()

