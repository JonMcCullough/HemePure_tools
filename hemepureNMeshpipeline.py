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

def transform_to_physical(pos, dx, shifts):
    return dx*(pos - shifts)

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
    xml += "  <visualisation>"
    xml += "     <centre units=\"m\" value=\"(0.0,0.0,0.0)\" />"
    xml += "     <orientation>"
    xml += "       <longitude units=\"deg\" value=\"45.0\" />"
    xml += "       <latitude units=\"deg\" value=\"45.0\" />"
    xml += "     </orientation>"
    xml += "     <display brightness=\"0.03\" zoom=\"1.0\" />"
    xml += "     <range>"
    xml += "       <maxvelocity units=\"m/s\" value=\"0.1\" />"
    xml += "       <maxstress units=\"Pa\" value=\"0.1\" />"
    xml += "     </range>"
    xml += "   </visualisation>"
    xml += "</hemelbsettings>\n";

    with open(hemexmlfname, "w") as outxml:
        outxml.write(xml)

## JM stuff ...
def write_dualMap(outN, inN):
    Aout = outletDict[outN]
    Vin = inletDict[inN]
    nI = len(Vin)
    nA = len(Aout)
    VinPaired = np.zeros(nI); VinPaired.fill(np.inf)
    dMinPaired = np.zeros(nI)
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
            dMinPaired[minID] = dmin
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
                dMinPaired[minID] = dmin
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
        dMinPaired[i] = dmin

    mapp = "Map Indices: Aout (index, radius) then Paired Vin(s) (index, distance, radius) \n"
    for o in range(nA):
        mapp += "(" + str(o) + "," + str(oRadDict[outN][o]) + ")"
        for i in range(nI):
            if VinPaired[i] == o: mapp += ", ( " + str(i) + "," + str(dMinPaired[i]) + "," + str(iRadDict[inN][i]) + ")"
        mapp += "\n"

    with open("mapAtoB.txt", "w") as outMap:
        outMap.write(mapp)

####
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

numMeshes = 2
inletDict={}
outletDict={}
iRadDict={}
oRadDict={}

for mesh in range(numMeshes):
    STLFNAME = sys.argv[1+mesh*5]
    STLUNITS = float(sys.argv[2+mesh*5])
    INLETS = [np.float_(iolet.split(",")) for iolet in (sys.argv[3+mesh*5]).split(";")]
    NUMINLETS = int(sys.argv[4+mesh*5])
    NUMOUTLETS = int(sys.argv[5+mesh*5])
    
    ROOTNAME = os.path.splitext(os.path.basename(STLFNAME))[0]

    
    
    print("Writing initial xml...")
    xmlfname = ROOTNAME + ".xml"

    inletpos0 = [np.array([0.0,0.0,0.0]) for i in range(NUMINLETS)]
    outletpos0 = [np.array([0.0,0.0,0.0]) for i in range(NUMOUTLETS)]
    write_voxelizer_xml(xmlfname, RESOLUTION, STLFNAME, inletpos0, outletpos0)

    # Run voxelizer but end early, dumping only the ioletpositions
    execute("mpirun -np " + str(NUMRANKS) + " " + VOXELIZERPATH + " " + xmlfname + "  ENDEARLY\n")
    execute("mv ioletpositions.txt ioletpositions_"+str(mesh)+".txt")

    iolet_list = []
    dx = None
    shifts = None
    with open("ioletpositions_"+str(mesh)+".txt", "r") as ioletpos:
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
    
    inletDict[mesh] = [transform_to_physical(i,dx,shifts) for i in inletposlist]
    outletDict[mesh] = [transform_to_physical(i,dx,shifts) for i in outletposlist]

    # Write the second version of the voxelizer's xml, in which the inlet and outlet positions are correctly identified and ordered
    write_voxelizer_xml(xmlfname, RESOLUTION, STLFNAME, inletposlist, outletposlist)

    # Run voxelizer to completion this time
    execute("mpirun -np " + str(NUMRANKS) + " " + VOXELIZERPATH + " " + xmlfname + "\n")
    execute("cat fluidAndLinks_*.dat > fluidAndLinks"+str(mesh)+".dat && rm fluidAndLinks_*.plb && rm fluidAndLinks_*.dat")
    execute("mv iolets_block_inputxml.txt iolets_block_inputxml_"+str(mesh)+".txt")

    with open("inlets_radius.txt","r") as ilets:
        iLETS = [line.rstrip('\n') for line in ilets]
        iLETS = [dx*float(i.split(',')[1]) for i in iLETS]
        ilets.close()

    iRadDict[mesh]=iLETS

    with open("outlets_radius.txt","r") as ilets:
        iLETS = [line.rstrip('\n') for line in ilets]
        iLETS = [dx*float(i.split(',')[1]) for i in iLETS]
        ilets.close()

    oRadDict[mesh]=iLETS

    execute("rm inlets_radius.txt")
    execute("rm outlets_radius.txt")

    # Get the inlets and outlets xml blocks (for the hemelb input xml) output by the voxelizer
    with open("iolets_block_inputxml_"+str(mesh)+".txt", "r") as ioletsblockfile:
        ioletsblocktxt = ioletsblockfile.read()

    # Write the hemelb input.xml file
    hemexmlfname = "input_"+str(mesh)+".xml"
    gmyfname = ROOTNAME + ".gmy"
    gmy_resolution = dx * STLUNITS
    write_heme_xml(hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt)

    # Convert the voxelizer output into a hemelb gmy file
    execute("bash " + MAKEGMYMPIPATH + " fluidAndLinks"+str(mesh)+".dat " + gmyfname  + " " + str(NUMRANKS) + " " + str(VX2GMY_CHUNKSIZE) + "\n")

# write mapping of outlets to inlets
print("Calculating and writing outlet/inlet mapping")
write_dualMap(0,1)

flow = "Velocity scale factors (*V_out): \n"
pressure = "dP scale factors (*P_out): \n"

with open("mapAtoB.txt") as theMap:
    next(theMap)
    for line in theMap:
        line = line.strip("\n").strip("(").strip(")").split("), (")
        outletData = [float(x) for x in line[0].split(",")]
        outletIDX = outletData[0]
        R02 = outletData[1]**2

        rInlets = np.zeros(len(line)-1)
        L_Inlets = np.zeros(len(line)-1)
        inletIDX = np.zeros(len(line)-1)
        for out in range(1,len(line)):
            inletData = [float(x) for x in line[out].split(",")]
            rInlets[out-1] = inletData[2]
            L_Inlets[out-1] = inletData[1]
            inletIDX[out-1] = inletData[0]
    
        DP = np.random.uniform(0.3, 0.7,len(rInlets)) ## Arbitrary calculation of pressure drop
        print("Pressure drop distribution = ", DP)
        R2 = rInlets**2
        
        print("r0 = ", R02)
        print("r = ", R2)
        print("L = ", L_Inlets)
        q = (R02/R2)*(DP/DP[0]*L_Inlets[0]/L_Inlets*R2/R2[0])*(np.sum(DP/DP[0]*L_Inlets[0]/L_Inlets*R2/R2[0]))**-1

        for i in range(len(rInlets)):
            flow += str(int(outletIDX)) + "," + str(int(inletIDX[i])) + "," + str(q[i]) + "\n"
            pressure += str(int(outletIDX)) + "," + str(int(inletIDX[i])) + "," + str(DP[i]) + "\n"

with open("velocityFactors.txt", "w") as factors:
    factors.write(flow)
factors.close()

with open("dPFactors.txt", "w") as factors:
    factors.write(pressure)
factors.close()




## Create the velocity weights file - WARNING: CURRENTLY ASSUMES ONLY 1 INLET (not easy to fix...)
#inletsfname = ROOTNAME + ".inlets"
#execute(GMY2INLETSPATH + " " + gmyfname + " " + inletsfname + "\n")
#execute("python " + INFLOWPROFILEBUILDERPATH + " " + inletsfname + "\n")
#execute("cp out.weights.txt INLET0_VELOCITY.txt.weights.txt\n")
