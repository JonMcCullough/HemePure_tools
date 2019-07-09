import os, sys, math

if len(sys.argv) != 10:
	sys.exit("python extract_specific_inlet.py input.inlets originx originy originz voxel_size posx posy posz radius\n")

infname = sys.argv[1]
originx = float(sys.argv[2])
originy = float(sys.argv[3])
originz = float(sys.argv[4])
voxel_size = float(sys.argv[5])
posx = float(sys.argv[6])
posy = float(sys.argv[7])
posz = float(sys.argv[8])
radius = float(sys.argv[9])

# calc what the real positions are in terms of lattice sites
lposx = (posx - originx)/voxel_size
lposy = (posy - originy)/voxel_size
lposz = (posz - originz)/voxel_size
lradius = radius/voxel_size

with open(infname, "r") as infile:
	for line in infile.readlines():
		sline = line.split()
		lx = float(sline[0])
		ly = float(sline[1])
		lz = float(sline[2])

		# Check if site is within lradius of specified position
		dx = lx - lposx
		dy = ly - lposy
		dz = lz - lposz
		dist = math.sqrt(dx*dx + dy*dy + dz*dz)

#		print lx, ly,lz, dist, "|", lposx, lposy, lposz, lradius

		if dist < lradius:
			print line,

