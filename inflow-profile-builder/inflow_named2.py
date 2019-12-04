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
  test_data = test_data.sort_values("iolet_number")
  test_data = test_data.reset_index(drop=True)
  adf = pd.DataFrame(test_data, columns=["x","y","z","iolet_number"])

  meshNum = sys.argv[2]
  
  #plot_plane.plot_xyz_plane_3views(tdf)
  border_points = test_data[test_data['site_type'] == 4] 
  bdf = pd.DataFrame(border_points, columns=["x","y","z","iolet_number"])

  inlet_points = test_data[test_data['site_type'] == 2] 
  idf = pd.DataFrame(inlet_points, columns=["x","y","z","iolet_number"])

  #import create_mesh_links
  #links = create_mesh_links.create_mesh_links(tdf, dist_limit=3.0)

  #border_points, border_point_ids = detect_borders.detect_border_points_in_plane(tdf_collapsed, links, threshold=14) #13 looked pretty good...
  #ax.scatter(border_points[:,0], border_points[:,1], border_points[:,2])
  #ax.scatter(bdf["x"], bdf["y"], bdf["z"])

  import calculate_distances
  distances = calculate_distances.calculate_distances_from_borders_per_inlet(adf, bdf, test_data['iolet_number'].max()+1)

  import distance_to_velocity_weights
  weights = distance_to_velocity_weights.distance_to_velocity_weights_per_inlet(distances)

  adf = adf.assign(nodeWeights=pd.Series(weights)) 
 
  for let in range(0,test_data['iolet_number'].max()+1):
      fig = plt.figure()
      ax = fig.add_subplot(121, projection='3d')
      ax.scatter(bdf[bdf['iolet_number'] == let]["x"], bdf[bdf['iolet_number'] == let]["y"], bdf[bdf['iolet_number'] == let]["z"])

      ax = fig.add_subplot(122, projection='3d')
      ax.scatter(adf[adf['iolet_number'] == let]["x"], adf[adf['iolet_number'] == let]["y"], adf[adf['iolet_number'] == let]["nodeWeights"])

      plt.title('This should be outline and weights of inlet plane' + str(let) + ' on Mesh' + meshNum)
      fig.savefig("InletImages/Mesh" + meshNum + "Inlet" + str(let) + "_results.png")
      plt.close()

      np.savetxt("out" + str(let) + ".weights.txt", adf[adf['iolet_number'] == let][['x','y','z','nodeWeights']].values, fmt='%d %d %d %f', delimiter=' ')

     

