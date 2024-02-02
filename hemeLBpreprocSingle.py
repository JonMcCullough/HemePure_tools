import os, sys
import numpy as np

VOXELIZERPATH = "./execs/voxelizer_MultiInput"
MAKEGMYMPIPATH = "./execs/make_gmy_MPI.sh"
VX2GMYPATH = "./execs/mpivx2gmy"

GMY2LETSPATH = "./execs/gmy2lets"
INFLOWPROFILEBUILDERPATH = "./execs/inflow-profile-builder/inflow_named3.py"

VX2GMY_CHUNKSIZE = 10000

def execute(command):
        print("Executing: " + command)
        r = os.system(command)
        if r != 0:
                sys.exit("Command failed.")

def transform_to_lattice(pos, dx, shifts):
    return pos/dx + shifts

def transform_to_physical(pos, dx, shifts):
    return dx*(pos - shifts)

def write_voxelizer_xml(xmlfname, dxREL, dxABS, STLFNAME, inletposlist, outletposlist):
    xml = '<?xml version="1.0" ?>\n<!-- the referenceDirection is used for the resolution -->\n<!-- see src/offLattice/triangularSurfaceMesh.hh -->\n<!-- 0 means x-direction, 1 means y-direction and 2 means z-direction -->\n<referenceDirection> 0 </referenceDirection>\n'
    xml += "<resolution> " + str(dxREL) + " </resolution>\n"
    xml += "<DX> " + str(dxABS) + " </DX>\n"
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

def write_heme_xml(tauValue, hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt, originShift):
    xml =  "<?xml version=\"1.0\"?>\n"
    xml += "<hemelbsettings version=\"3\">\n"
    xml += "  <simulation>\n"
    xml += "    <step_length units=\"s\" value=\"" + str((tauValue - 0.5)*gmy_resolution*gmy_resolution/(12e-6)) + "\"/>\n"
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
    xml += "  <properties>\n"
    xml += "   <propertyoutput file=\"inlet.dat\" period=\"100\">\n"
    xml += "     <geometry type=\"inlet\" />\n"
    xml += "     <field type=\"velocity\" />\n"
    xml += "     <field type=\"pressure\" />\n"
    xml += "   </propertyoutput>\n"
    xml += "   <propertyoutput file=\"outlet.dat\" period=\"100\">\n"
    xml += "     <geometry type=\"outlet\" />\n"
    xml += "     <field type=\"velocity\" />\n"
    xml += "     <field type=\"pressure\" />\n"
    xml += "   </propertyoutput>\n"
    xml += "   <propertyoutput file=\"whole.dat\" period=\"100\">\n"
    xml += "     <geometry type=\"whole\" />\n"
    xml += "     <field type=\"velocity\" />\n"
    xml += "     <field type=\"pressure\" />\n"
    xml += "   </propertyoutput>\n"
    xml += "  </properties>\n"
    xml += "</hemelbsettings>\n";

    with open(hemexmlfname, "w") as outxml:
        outxml.write(xml)


if len(sys.argv) != 10:
    sys.exit("Usage: python3 hemeLBPreProc.py STLFNAME STLUNITS(e.g 1e-3 for mm) INLETPOSITIONS(X1,Y1,Z1;X2,Y2,Z2;..., (in quotes)) NUMINLETS NUMOUTLETS DXreq NUMRANKS RANKSPERNODE tauDesired")

NUMRANKS = int(sys.argv[-3])
RANKSPERNODE = int(sys.argv[-2])
CORESPERRANK=64/RANKSPERNODE
tauDes = float(sys.argv[-1])

STLFNAME = sys.argv[1]
STLUNITS = float(sys.argv[2])
#INLETS = [np.float_(iolet.split(",")) for iolet in (sys.argv[3]).split(";")]
with open(sys.argv[3]) as inletList:
    INLETS = [np.float_(iolet.split(",")) for iolet in inletList.readline().split(";")]
NUMINLETS = int(sys.argv[4])
NUMOUTLETS = int(sys.argv[5])
DXreq = np.float64(sys.argv[6])

ROOTNAME = os.path.splitext(os.path.basename(STLFNAME))[0]

print("STLFNAME = ", STLFNAME)
print("STLUNITS = ", STLUNITS)
print("INLETS = ", INLETS)
print("NUMINLETS = ", NUMINLETS)
print("NUMOUTLETS = ", NUMOUTLETS)
print("DXreq = ", DXreq)

print("ROOTNAME = ", ROOTNAME)

print("Writing initial xml...")
xmlfname = ROOTNAME + ".xml"

inletpos0 = [np.array([0.0,0.0,0.0]) for i in range(NUMINLETS)]
outletpos0 = [np.array([0.0,0.0,0.0]) for i in range(NUMOUTLETS)]
write_voxelizer_xml(xmlfname, DXreq/STLUNITS, DXreq, STLFNAME, inletpos0, outletpos0)

# Run voxelizer but end early, dumping only the ioletpositions
execute("srun -n " + str(int(NUMRANKS)) + " -N  " + str(int(NUMRANKS/RANKSPERNODE)) + " -c "+str(int(CORESPERRANK))+"  --unbuffered "  + VOXELIZERPATH + " " + xmlfname + "  ENDEARLY\n")

iolet_list = []
dx = None
shifts = None
with open("ioletpositions.txt", "r") as ioletpos:
    lines = ioletpos.readlines()
    
    # Get the real resolution
    if not lines[0].startswith('DX:'):
        sys.exit("ioletpositions.txt output from voxelizer does not have DX: line where expected (first line)")
    dx = float(lines[0].split()[1])
    print("dx RELATIVE = ", dx)
    dx = DXreq #dx*STLUNITS
    print("dx ABSOLUTE = ", dx)

    # Get the shift the voxelizer is applying to the STL
    if not lines[1].startswith('SHIFTS:'):
        sys.exit("ioletpositions.txt output from voxelizer does not have DX: line where expected (first line)")
    shifts = np.array([float(i) for i in lines[1].split()[1:]])
    print("shifts = ", shifts)

    ioletPosPU = lines[0]
    ioletPosPU += lines[1]
    # Get the iolet positions
    for line in lines[2:]:
        ioletposFull = [float(i) for i in line.split()]
        iolet_list.append(np.array(ioletposFull[0:3]))
        ioletpos = transform_to_physical(np.array(ioletposFull[0:3]),dx,shifts)
        ioletPosPU += str(ioletpos[0]) + " " + str(ioletpos[1]) + " " + str(ioletpos[2]) + "\n"
        #ioletPosPU += str(ioletpos[0]) + " " + str(ioletpos[1]) + " " + str(ioletpos[2]) + " " + str(ioletposFull[3]*dx)+ "\n"

with open("ioletpositions_PU.txt", "w") as outxml:
    outxml.write(ioletPosPU)
print("Written iolets file in PU, may be handy for inlets text file")

# Work out the inlet positions (provided to this script) in lattice units
INLETS_LATTICE = [transform_to_lattice(inletpos, dx, shifts) for inletpos in INLETS]

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

#raise SystemExit(0)

# Write the second version of the voxelizer's xml, in which the inlet and outlet positions are correctly identified and ordered
write_voxelizer_xml(xmlfname, DXreq/STLUNITS, DXreq, STLFNAME, inletposlist, outletposlist)

# Run voxelizer to completion this time
execute("srun -n " + str(int(NUMRANKS)) + " -N " + str(int(NUMRANKS/RANKSPERNODE)) + " -c "+str(int(CORESPERRANK))+"  --unbuffered " + VOXELIZERPATH + " " + xmlfname + "\n")

# Catting many files should no longer be required. The voxelizer should dump just one large file.
execute("rm -f fluidsAndLinks.plb")

# Write the hemelb input.xml files - with different BCs
gmyfname = ROOTNAME + ".gmy"
gmy_resolution = dx #* STLUNITS

#BCExtensions=["PP","VP","VfP","VfWKf"]
BCExtensions=["PP"]

for t in range(len(BCExtensions)):
    hemexmlfname = "input_"+BCExtensions[t]+".xml"

    # Get the inlets and outlets xml blocks (for the hemelb input xml) output by the voxelizer
    with open("iolets_block_inputxml_"+BCExtensions[t]+".txt", "r") as ioletsblockfile:
        ioletsblocktxt = ioletsblockfile.read()
    write_heme_xml(tauDes, hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt, dx*shifts) #JM had - shiftMaster here

# Convert the voxelizer output into a hemelb gmy file
execute("srun -n " + str(int(NUMRANKS)) + " -N  " + str(int(NUMRANKS/RANKSPERNODE)) + " -c "+str(int(CORESPERRANK))+"  --unbuffered "  + VX2GMYPATH + " ./fluidsAndLinks.dat __make_gmy_MPI_test.gmy " + str(VX2GMY_CHUNKSIZE) + "\n")
execute("cat __make_gmy_MPI_test.gmy __make_gmy_MPI_test.gmy_blockdata_temp_ > "+ gmyfname+"\n")

## Create the velocity weights file 
if not os.path.exists('InletImages'):
    os.mkdir('InletImages')
inletsfname = ROOTNAME + ".inlets"
execute(GMY2LETSPATH + " " + gmyfname + " " + inletsfname + " INLET \n")
execute("python3 " + INFLOWPROFILEBUILDERPATH + " " + inletsfname + " 0 INLET \n")

for ilet in range(0,NUMINLETS): 
    execute("cp out" + str(ilet) + ".txt.weights.txt INLET" + str(ilet) + "_VELOCITY.txt.weights.txt\n") 

## Create the windkessel weights file 
if not os.path.exists('OutletImages'):
    os.mkdir('OutletImages')
outletsfname = ROOTNAME + ".outlets"
execute(GMY2LETSPATH + " " + gmyfname + " " + outletsfname + " OUTLET \n")
execute("python3 " + INFLOWPROFILEBUILDERPATH + " " + outletsfname + " 0 OUTLET \n")

for ilet in range(0,NUMOUTLETS): 
    execute("cp out" + str(ilet) + ".txt.weights.txt OUTLET" + str(ilet) + "_WK.txt.weights.txt\n") 

