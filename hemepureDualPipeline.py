import os, sys
import numpy as np

VOXELIZERPATH = "/cs/heme/HemePure_JM/HemePure_tools/voxelizer/source/voxelizer "
MAKEGMYMPIPATH = "/cs/heme/HemePure_JM/HemePure_tools/vx2gmy/make_gmy_MPI.sh"
GMY2INLETSPATH = "~/gmyTools/gmy2inlets/gmy2inlets"
INFLOWPROFILEBUILDERPATH = "~/inflow-profile-builder/inflow.py"

NUMRANKS = 2
VX2GMY_CHUNKSIZE = 2000

def execute(command):
        print("Executing: " + command)
        r = os.system(command)
        if r != 0:
                sys.exit("Command failed.")

def transform_to_lattice(pos, dx, shifts):
    return pos/dx + shifts

def write_voxelizer_xml(xmlfname, RESOLUTION, STLFNAME, inletposlist, outletposlist):
    xml = '<?xml version="1.0" ?>\n<!-- the referenceDirection is used for the resolution -->\n<!-- see src/offLattice/triangularSurfaceMesh.hh -->\n<!-- 0 means x-direction, 1 means y-direction and 2 means z-direction -->\n<referenceDirection> 0 </referenceDirection>\n'
    xml += "<resolution> " + str(RESOLUTION) + " </resolution>\n"
    xml += "<!-- *.stl containing geometry -->\n"
    xml += "<stl> " + STLFNAME + " </stl>\n"
    xml += "<!-- analysis points for identification of iolets -->\n<!-- first <num_Ilets> points identify inlets -->\n<!-- last <num_Olets> points identify outlets -->\n<analysisPoints> <!-- lattice units -->\n"
    xml += "<numIlets> " + str(len(inletposlist)) + " </numIlets>\n"
    xml += "<numOlets> " + str(len(outletposlist)) + " </numOlets>\n"

    iolet = 1
    for pos in inletposlist:
        xml += '<point id="'+str(iolet)+'"> '
        xml += str(pos[0]) + ' ' + str(pos[1]) + ' ' + str(pos[2])
        xml += ' </point>\n'
        iolet += 1
    for pos in outletposlist:
        xml += '<point id="'+str(iolet)+'"> '
        xml += str(pos[0]) + ' ' + str(pos[1]) + ' ' + str(pos[2])
        xml += ' </point>\n'
        iolet += 1
    xml += "</analysisPoints>\n"

    with open(xmlfname, "w") as outxml:
        outxml.write(xml)

def write_heme_xml(hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt):
    xml =  "<?xml version=\"1.0\"?>\n"
    xml += "<hemelbsettings version=\"3\">\n"
    xml += "  <simulation>\n"
    xml += "    <step_length units=\"s\" value=\"CHANGE\"/>\n"
    xml += "    <steps units=\"lattice\" value=\"CHANGE\"/>\n"
    xml += "    <stresstype value=\"1\"/>\n"
    xml += "    <voxel_size units=\"m\" value=\"" + str(gmy_resolution) + "\"/>\n"
    xml += "    <origin units=\"m\" value=\"(0.0,0.0,0.0)\"/>\n"
    xml += "  </simulation>\n"
    xml += " <geometry>\n"
    xml += "    <datafile path=\"" + gmyfname + "\"/>\n"
    xml += "  </geometry>\n"
    xml += "  <initialconditions>\n"
    xml += "    <pressure>\n"
    xml += "      <uniform units=\"mmHg\" value=\"0.0\"/>\n"
    xml += "    </pressure>\n"
    xml += "  </initialconditions>\n"
    xml += "  <monitoring>\n"
    xml += "    <incompressibility/>\n"
    xml += "  </monitoring>\n\n"
    xml += ioletsblocktxt + "\n"
    xml += "  <visualisation>\n"
    xml += "     <centre units=\"m\" value=\"(0.0,0.0,0.0)\" />\n"
    xml += "     <orientation>\n"
    xml += "       <longitude units=\"deg\" value=\"45.0\" />\n"
    xml += "       <latitude units=\"deg\" value=\"45.0\" />\n"
    xml += "     </orientation>\n"
    xml += "     <display brightness=\"0.03\" zoom=\"1.0\" />\n"
    xml += "     <range>\n"
    xml += "       <maxvelocity units=\"m/s\" value=\"0.1\" />\n"
    xml += "       <maxstress units=\"Pa\" value=\"0.1\" />\n"
    xml += "     </range>\n"
    xml += "   </visualisation>\n"
    xml += "</hemelbsettings>\n";

    with open(hemexmlfname, "w") as outxml:
        outxml.write(xml)


def write_dualMap(Aout, Vin):
    nI = len(Vin)
    nA = len(Aout)
    VinPaired = np.zeros(nI); VinPaired.fill(np.inf)
    
    k=[]
    
    for o in range(nA):
        dmin = 5
        minID = -1
        for i in range(nI):
            dnow = np.linalg.norm(Aout[o] - Vin[i])
            if dnow < dmin:
                dmin = dnow
                minID = i
              
        if not minID in k:
            VinPaired[minID] = o
            k.append(minID)
        else:
                sub = np.delete(np.arange(nI),k)
                
                dmin = 5
                minID = -1
                for i in sub:
                    dnow = np.linalg.norm(Aout[o] - Vin[i])
                    if dnow < dmin:
                        dmin = dnow
                        minID = i
                
                VinPaired[minID] = o
                k.append(minID)

    unVI = np.delete(np.arange(nI),k)

    for i in unVI:
        dmin = float("+Inf")
        minID = -1
        for o in range(nA):
            dnow = np.linalg.norm(Aout[o] - Vin[i])
            if dnow < dmin:
                dmin = dnow
                minID = o
        VinPaired[i] = minID
    
    mapp = "Map Indices: Aout, Paired Vin(s) \n"
    for o in range(nA):
        mapp += str(o)
        for i in range(nI):
            if VinPaired[i] == o: mapp += "," + str(i)
        mapp += "\n"
    
    with open("mapAtoV.txt", "w") as outMap:
        outMap.write(mapp)




if len(sys.argv) != 12:
    sys.exit("Usage: python3 hemepureDualpipeline.py  STLFNAME_A STLUNITS_A(e.g 1e-3 for mm) INLETPOSITIONS_A(X1,Y1,Z1;X2,Y2,Z2;..., (in quotes)) NUMINLETS_A NUMOUTLETS_A STLFNAME_B STLUNITS_B(e.g 1e-3 for mm) INLETPOSITIONS_B(X1,Y1,Z1;X2,Y2,Z2;... (in quotes)) NUMINLETS_B NUMOUTLETS_B RESOLUTION(stupid palabos units e.g. 150)")

STLFNAME_A = sys.argv[1]
STLUNITS_A = float(sys.argv[2])
INLETS_A = [np.float_(iolet.split(",")) for iolet in (sys.argv[3]).split(";")]
NUMINLETS_A = int(sys.argv[4])
NUMOUTLETS_A = int(sys.argv[5])

STLFNAME_B = sys.argv[6]
STLUNITS_B = float(sys.argv[7])
INLETS_B = [np.float_(iolet.split(",")) for iolet in (sys.argv[8]).split(";")]
NUMINLETS_B = int(sys.argv[9])
NUMOUTLETS_B = int(sys.argv[10])

RESOLUTION = int(sys.argv[11])

ROOTNAME_A = os.path.splitext(os.path.basename(STLFNAME_A))[0]
ROOTNAME_B = os.path.splitext(os.path.basename(STLFNAME_B))[0]

print("Arterial Tree Values")
print("STLFNAME = ", STLFNAME_A)
print("STLUNITS = ", STLUNITS_A)
print("INLETS = ", INLETS_A)
print("NUMINLETS = ", NUMINLETS_A)
print("NUMOUTLETS = ", NUMOUTLETS_A)
print("RESOLUTION = ", RESOLUTION)

print("ROOTNAME = ", ROOTNAME_A)


print("Venous Tree Values")
print("STLFNAME = ", STLFNAME_B)
print("STLUNITS = ", STLUNITS_B)
print("INLETS = ", INLETS_B)
print("NUMINLETS = ", NUMINLETS_B)
print("NUMOUTLETS = ", NUMOUTLETS_B)
print("RESOLUTION = ", RESOLUTION)

print("ROOTNAME = ", ROOTNAME_B)

print("Writing initial xml...")
xmlfname_A = ROOTNAME_A + ".xml"
xmlfname_B = ROOTNAME_B + ".xml"

inletpos0 = [np.array([0.0,0.0,0.0]) for i in range(NUMINLETS_A)]
outletpos0 = [np.array([0.0,0.0,0.0]) for i in range(NUMOUTLETS_A)]
write_voxelizer_xml(xmlfname_A, RESOLUTION, STLFNAME_A, inletpos0, outletpos0)

inletpos1 = [np.array([0.0,0.0,0.0]) for i in range(NUMINLETS_B)]
outletpos1 = [np.array([0.0,0.0,0.0]) for i in range(NUMOUTLETS_B)]
write_voxelizer_xml(xmlfname_B, RESOLUTION, STLFNAME_B, inletpos1, outletpos1)

# Run voxelizer but end early, dumping only the ioletpositions
execute("mpirun -np " + str(NUMRANKS) + " " + VOXELIZERPATH + " " + xmlfname_A + "  ENDEARLY\n")
execute("mpirun -np " + str(NUMRANKS) + " " + VOXELIZERPATH + " " + xmlfname_B + "  ENDEARLY\n")

### MESH A
print("\n IOlets for Mesh A...")
iolet_list = []
dx = None
shifts = None
with open("ioletpositions.txt", "r") as ioletpos:
    lines = ioletpos.readlines()

    # Get the real resolution
    if not lines[0].startswith('DX:'):
        sys.exit("ioletpositions.txt output from voxelizer does not have DX: line where expected (first line)")
    dx = float(lines[0].split()[1])
    print("dx = ", dx)

    # Get the shift the voxelizer is applying to the STL
    if not lines[1].startswith('SHIFTS:'):
        sys.exit("ioletpositions.txt output from voxelizer does not have DX: line where expected (first line)")
    shifts = np.array([float(i) for i in lines[1].split()[1:]])
    print("shifts = ", shifts)

    # Get the iolet positions
    for line in lines[2:]:
        ioletpos = [float(i) for i in line.split()]
        iolet_list.append(np.array(ioletpos))

# Work out the inlet positions (provided to this script) in lattice units
INLETS_LATTICE = [transform_to_lattice(inletpos, dx, shifts) for inletpos in INLETS_A]

# Identify the closest iolets to the iolet positions passed to this script
inlets_list = []
for inletpos in INLETS_LATTICE:
    min_dist = float("+Inf")
    favoured_ioindex = -1
    for ioindex, ioletpos in enumerate(iolet_list):
        dist = np.linalg.norm(ioletpos - inletpos)
        if dist < min_dist:
            min_dist = dist
            favoured_ioindex = ioindex
    # Check that the closest inlet is not already in the list (for a different inlet)
    # This would suggest that the user has entered wrong positions (or two openings are
    # ridiculously close to each other)
    if favoured_ioindex in inlets_list:
        sys.exit("inletpos " + str(inletpos) + " corresponds to more than one 'nearest' opening")
    inlets_list.append(favoured_ioindex)

print("Identified inlet(s) by index:")
print(inlets_list)

inletposlist = []
outletposlist = []
for ioindex, ioletpos in enumerate(iolet_list):

    # If index is not in list of inlet indices then it's an outlet
    if ioindex not in inlets_list:
        outletposlist.append(ioletpos)
    else:
        inletposlist.append(ioletpos)

outletsA = outletposlist
# Write the second version of the voxelizer's xml, in which the inlet and outlet positions are correctly identified and ordered
write_voxelizer_xml(xmlfname_A, RESOLUTION, STLFNAME_A, inletposlist, outletposlist)


### MESH B
print("\n IOlets for Mesh B...")
iolet_list = []
dx = None
shifts = None
with open("ioletpositions.txt", "r") as ioletpos:
    lines = ioletpos.readlines()

    # Get the real resolution
    if not lines[0].startswith('DX:'):
        sys.exit("ioletpositions.txt output from voxelizer does not have DX: line where expected (first line)")
    dx = float(lines[0].split()[1])
    print("dx = ", dx)

    # Get the shift the voxelizer is applying to the STL
    if not lines[1].startswith('SHIFTS:'):
        sys.exit("ioletpositions.txt output from voxelizer does not have DX: line where expected (first line)")
    shifts = np.array([float(i) for i in lines[1].split()[1:]])
    print("shifts = ", shifts)

    # Get the iolet positions
    for line in lines[2:]:
        ioletpos = [float(i) for i in line.split()]
        iolet_list.append(np.array(ioletpos))

# Work out the inlet positions (provided to this script) in lattice units
INLETS_LATTICE = [transform_to_lattice(inletpos, dx, shifts) for inletpos in INLETS_B]

# Identify the closest iolets to the iolet positions passed to this script
inlets_list = []
for inletpos in INLETS_LATTICE:
    min_dist = float("+Inf")
    favoured_ioindex = -1
    for ioindex, ioletpos in enumerate(iolet_list):
        dist = np.linalg.norm(ioletpos - inletpos)
        print("checking")
        print(ioindex,dist)
        if dist < min_dist:
            min_dist = dist
            favoured_ioindex = ioindex
    # Check that the closest inlet is not already in the list (for a different inlet)
    # This would suggest that the user has entered wrong positions (or two openings are
    # ridiculously close to each other)
    if favoured_ioindex in inlets_list:
        sys.exit("inletpos " + str(inletpos) + " corresponds to more than one 'nearest' opening")
    inlets_list.append(favoured_ioindex)

print("Identified inlet(s) by index:")
print(inlets_list)

inletposlist = []
outletposlist = []
for ioindex, ioletpos in enumerate(iolet_list):

    # If index is not in list of inlet indices then it's an outlet
    if ioindex not in inlets_list:
        outletposlist.append(ioletpos)
    else:
        inletposlist.append(ioletpos)

inletsB = inletposlist;
# Write the second version of the voxelizer's xml, in which the inlet and outlet positions are correctly identified and ordered
write_voxelizer_xml(xmlfname_B, RESOLUTION, STLFNAME_B, inletposlist, outletposlist)

print("outlets A", outletsA)
print("inlets B", inletsB)

# Run voxelizer to completion this time
execute("mpirun -np " + str(NUMRANKS) + " " + VOXELIZERPATH + " " + xmlfname_A + "\n")
execute("cat fluidAndLinks_*.dat > fluidAndLinksA.dat && rm fluidAndLinks_*.plb && rm fluidAndLinks_*.dat")
execute("mv iolets_block_inputxml.txt iolets_block_inputxml_A.txt")

execute("mpirun -np " + str(NUMRANKS) + " " + VOXELIZERPATH + " " + xmlfname_B + "\n")
execute("cat fluidAndLinks_*.dat > fluidAndLinksB.dat && rm fluidAndLinks_*.plb && rm fluidAndLinks_*.dat")
execute("mv iolets_block_inputxml.txt iolets_block_inputxml_B.txt")

### MESH A
# Get the inlets and outlets xml blocks (for the hemelb input xml) output by the voxelizer
with open("iolets_block_inputxml_A.txt", "r") as ioletsblockfile:
    ioletsblocktxt = ioletsblockfile.read()

# Write the hemelb input.xml file
hemexmlfname_A = "input_A.xml"
gmyfname_A = ROOTNAME_A + ".gmy"
gmy_resolution = dx * STLUNITS_A
write_heme_xml(hemexmlfname_A, gmyfname_A, gmy_resolution, ioletsblocktxt)

# Convert the voxelizer output into a hemelb gmy file
execute("bash " + MAKEGMYMPIPATH + " fluidAndLinksA.dat " + gmyfname_A  + " " + str(NUMRANKS) + " " + str(VX2GMY_CHUNKSIZE) + "\n")


### MESH B
# Get the inlets and outlets xml blocks (for the hemelb input xml) output by the voxelizer
with open("iolets_block_inputxml_B.txt", "r") as ioletsblockfile:
    ioletsblocktxt = ioletsblockfile.read()

# Write the hemelb input.xml file
hemexmlfname_B = "input_B.xml"
gmyfname_B = ROOTNAME_B + ".gmy"
gmy_resolution = dx * STLUNITS_B
write_heme_xml(hemexmlfname_B, gmyfname_B, gmy_resolution, ioletsblocktxt)

# Convert the voxelizer output into a hemelb gmy file
execute("bash " + MAKEGMYMPIPATH + " fluidAndLinksB.dat " + gmyfname_B  + " " + str(NUMRANKS) + " " + str(VX2GMY_CHUNKSIZE) + "\n")

# write mapping of outlets to inlets
print("Calculating and writing outlet/inlet mapping")
write_dualMap(outletsA,inletsB)

## Create the velocity weights file - WARNING: CURRENTLY ASSUMES ONLY 1 INLET (not easy to fix...)
#inletsfname = ROOTNAME + ".inlets"
#execute(GMY2INLETSPATH + " " + gmyfname + " " + inletsfname + "\n")
#execute("python " + INFLOWPROFILEBUILDERPATH + " " + inletsfname + "\n")
#execute("cp out.weights.txt INLET0_VELOCITY.txt.weights.txt\n")
