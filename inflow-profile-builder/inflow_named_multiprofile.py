import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import fit_plane
import plot_plane
import transform_points
import detect_borders
import os, sys


if __name__ == "__main__":

  #test_data = pd.read_csv("test.csv", delim_whitespace=True)
  test_data = pd.read_csv(sys.argv[1], delim_whitespace=True)
  test_data.columns = ["x","y","z","site_type","iolet_number"]
  test_data = test_data.sort_values("iolet_number")
  test_data = test_data.reset_index(drop=True)
  adf = pd.DataFrame(test_data, columns=["x","y","z","iolet_number"])

  meshNum = sys.argv[2]
  
  letType = str(sys.argv[3])

  plotBool = False 
  if len(sys.argv)==5:
      plotBool = bool(sys.argv[4])

  if (letType == 'INLET'):
      #plot_plane.plot_xyz_plane_3views(tdf)
      border_points = test_data[test_data['site_type'] == 4] 
      bdf = pd.DataFrame(border_points, columns=["x","y","z","iolet_number"])

      inlet_points = test_data[test_data['site_type'] == 2] 
      idf = pd.DataFrame(inlet_points, columns=["x","y","z","iolet_number"])
  elif (letType == 'OUTLET'):
      #plot_plane.plot_xyz_plane_3views(tdf)
      border_points = test_data[test_data['site_type'] == 5] 
      bdf = pd.DataFrame(border_points, columns=["x","y","z","iolet_number"])

      inlet_points = test_data[test_data['site_type'] == 3] 
      idf = pd.DataFrame(inlet_points, columns=["x","y","z","iolet_number"])
  else:
      sys.exit('INVALID IOlet type: INLET or OUTLET')


  #import create_mesh_links
  #links = create_mesh_links.create_mesh_links(tdf, dist_limit=3.0)

  #border_points, border_point_ids = detect_borders.detect_border_points_in_plane(tdf_collapsed, links, threshold=14) #13 looked pretty good...
  #ax.scatter(border_points[:,0], border_points[:,1], border_points[:,2])
  #ax.scatter(bdf["x"], bdf["y"], bdf["z"])

  import calculate_distances
  distances = calculate_distances.calculate_distances_from_borders_per_inlet(adf, bdf, test_data['iolet_number'].max()+1)

  import distance_to_velocity_weights
  weights_poiseuille = distance_to_velocity_weights.distance_to_velocity_weights_per_inlet(distances)
  weights_plug = distance_to_velocity_weights.distance_to_plugvelocity_weights_per_inlet(distances)
  weights_stenosis = distance_to_velocity_weights.distance_to_stenosisvelocity_weights_per_inlet(distances)

  weightsets = {'Poiseuille': weights_poiseuille,
          'Plug': weights_plug,
          'Stenosis':weights_stenosis}

  for name,wts in weightsets.items():
      adf = adf.assign(nodeWeights=pd.Series(wts)) 
 
      for let in range(0,test_data['iolet_number'].max()+1):
           np.savetxt("out" + str(let) + "_" + str(name) + ".txt.weights.txt", adf[adf['iolet_number'] == let][['x','y','z','nodeWeights']].values, fmt='%d %d %d %f', delimiter=' ') 
      
  #Generate plots of each inlet profile if desired
  if plotBool:
      if letType=='INLET' and not os.path.exists('InletImages'):
          os.mkdir('InletImages')
      if letType=='OUTLET' and not os.path.exists('OutletImages'):
          os.mkdir('OutletImages')
      
      for let in range(0,test_data['iolet_number'].max()+1):
          temp = bdf[bdf['iolet_number'] == let]
          localPoints = temp[['x','y','z']].values
      
          # subtract out the centroid   
          localPoints = localPoints - localPoints.mean(axis=0) 
    
          # singular value decomposition
          u,sig,v = np.linalg.svd(localPoints)
          normal = v[2]
          print(v[2])    
          azim = np.arccos(np.dot(np.array([normal[0],normal[1],0]),np.array([1,0,0])) / (np.linalg.norm(np.array([normal[0],normal[1],0])) * np.linalg.norm(np.array([1,0,0]))))
          elev = np.pi/2 - np.arccos(np.dot(normal,np.array([0,0,1])) / (np.linalg.norm(normal) * np.linalg.norm(np.array([0,0,1]))))
          
          if np.isnan(azim):
              azim = 0

          print(np.degrees(azim), np.degrees(elev))

          fig = plt.figure()
          ax = fig.add_subplot(121, projection='3d')
          ax.scatter(bdf[bdf['iolet_number'] == let]["x"], bdf[bdf['iolet_number'] == let]["y"], bdf[bdf['iolet_number'] == let]["z"])
          ax.azim = np.degrees(azim)
          ax.elev = np.degrees(elev)

          ax = fig.add_subplot(122, projection='3d')
          ax.scatter(adf[adf['iolet_number'] == let]["x"], adf[adf['iolet_number'] == let]["y"], adf[adf['iolet_number'] == let]["nodeWeights"])

          if letType == 'INLET':
              plt.title('This should be outline and weights of inlet plane' + str(let) + ' on Mesh' + meshNum)
              fig.savefig("InletImages/Mesh" + meshNum + "Inlet" + str(let) + "_results.png")
          else:
              plt.title('This should be outline and weights of outlet plane' + str(let) + ' on Mesh' + meshNum)
              fig.savefig("InletImages/Mesh" + meshNum + "Outlet" + str(let) + "_results.png")

          plt.close()



     

