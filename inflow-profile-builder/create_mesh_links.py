import numpy as np
import pandas as pd
import sys

def create_mesh_links(tdf,dist_limit=1.0):
  """ Take a point cloud, and collapse it to the plane specified.
  """

  links = []

  px = []
  py = []
  pz = []

  num_sites = len(tdf["x"]) 
  interval = int(np.sqrt(num_sites))

  #convert to matrix for better performance
  xs = pd.DataFrame(tdf, columns=['px', 'py', 'pz']).to_numpy()

  for i in range(0,num_sites): #JM was xrange
    links_i = []
    if i%interval == 0:
      print i

    xsi1 = 1 + xs[i][0]
    yi = xs[i][1]
    zi = xs[i][2] 
    # O(N^2) region    
    for j in range(0,num_sites): #JM was xrange
      r = xsi1 - xs[j][0] 

      #replaced because it's 50% slower
      #if np.abs(xs[i][0] - xs[j][0]) < dist_limit:

      if r>=0:
        if r<=2: 
          if np.abs(yi - xs[j][1]) < dist_limit: 
            if np.abs(zi - xs[j][2]) < dist_limit: 

              # 1-liner to calculate distance between two points.
              dist = np.linalg.norm(xs[i]-xs[j])

              if dist < dist_limit:
                if i != j:
                  links_i.append(j)
    links.append(links_i)

  print "diagnostic: printing all lattice sites which have <2 links to neighbours"
  for i in range(0, len(links)): #JM was xrange
    if(len(links[i]) < 2):
      print i, links[i]
      
  return links

