import os, sys
import numpy as np

VOXELIZERPATH = "/cs/heme/HemePure_JM/HemePure_tools/voxelizer/source/voxelizer2"
MAKEGMYMPIPATH = "/cs/heme/HemePure_JM/HemePure_tools/vx2gmy/make_gmy_MPI.sh"
GMY2INLETSPATH = "~/gmyTools/gmy2inlets/gmy2inlets"
INFLOWPROFILEBUILDERPATH = "~/inflow-profile-builder/inflow.py"

# NUMRANKS = 6
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

def write_voxelizer_xml(xmlfname, DXreq, STLFNAME, inletposlist, outletposlist):
    xml = '<?xml version="1.0" ?>\n<!-- the referenceDirection is used for the resolution -->\n<!-- see src/offLattice/triangularSurfaceMesh.hh -->\n<!-- 0 means x-direction, 1 means y-direction and 2 means z-direction -->\n<referenceDirection> 0 </referenceDirection>\n'
    xml += "<DX> " + str(DXreq) + " </DX>\n"
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

def write_heme_xml(meshID, tauValue, hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt, originShift):
    xml =  "<?xml version=\"1.0\"?>\n"
    xml += "<hemelbsettings version=\"3\">\n"
    xml += "<!-- Mesh Number = " + str(meshID) + " -->\n"
    xml += "  <simulation>\n"
    xml += "    <step_length units=\"s\" value= " + str((tauValue - 0.5)*gmy_resolution*gmy_resolution/(12e-6)) + "/>\n"
    xml += "    <steps units=\"lattice\" value=\"CHANGE\"/>\n"
    xml += "    <stresstype value=\"1\"/>\n"
    xml += "    <voxel_size units=\"m\" value=\"" + str(gmy_resolution) + "\"/>\n"
    xml += "    <origin units=\"m\" value=\"(0.0,0.0,0.0)\"/>\n"
    xml += "  </simulation>\n"
    xml += " <geometry>\n"
    xml += "    <datafile path=\"" + gmyfname + "\"/>\n"
    xml += "    <mapping path=\"scalingMapXtoY.txt\" value=\"" + str(meshID) + "\"/>\n"
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
    xml += "   <propertyoutput file=\"coupled\" period=\"100\">\n"
    xml += "    <geometry type=\"XXXlet\" />\n"
    xml += "     <field type=\"coupled\" />\n"
    xml += "   </propertyoutput>\n"
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
    xml += "  </properties>\n"
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

    with open("mapMesh" + str(outN) + "to" + str(inN) + ".txt", "w") as outMap:
        outMap.write(mapp)


def write_scalingMapFile(outN, inN):

    #scalefile = "Scale file: outletMeshID outletIDX inletMeshIDX (inletIDX velocityScaleFactor pressureScaleFactor)_repeatedForAllInlets: \n"
    scalefile=""
    
    with open("mapMesh" + str(outN) + "to" + str(inN) + ".txt") as theMap:
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
        
            R2 = rInlets**2
            
            print("r0 = ", R02)
            print("r = ", R2)
            print("L = ", L_Inlets)
            
            if (len(rInlets) == 1 and L_Inlets[0] == 0.0):
                DP = np.zeros(len(rInlets))
                q = np.ones(len(rInlets))
                print("Special case - direct mapping (1:1) between meshes -> DP = 0, q = 1")
                print("Pressure drop distribution = ", DP)
                print("Scale factors for velocity = ", q)

            else:
                DP = np.random.uniform(0.3, 0.7,len(rInlets)) ## Arbitrary calculation of pressure drop
                q = (R02/R2)*(DP/DP[0]*L_Inlets[0]/L_Inlets*R2/R2[0])*(np.sum(DP/DP[0]*L_Inlets[0]/L_Inlets*R2/R2[0]))**-1 
                print("Pressure drop distribution = ", DP)
                print("Scale factors for velocity = ", q)


            scalefile += str(int(outN)) + "," + str(int(outletIDX)) + "," + str(int(inN))

            for i in range(len(rInlets)):
                scalefile += "," + str(int(inletIDX[i])) + "," + str(q[i]) + "," + str(DP[i]) 
            scalefile += "\n"
    theMap.close()

    with open("scalingMap" + str(outN) + "to" + str(inN) + ".txt", "w") as factors:
        factors.write(scalefile)
    factors.close()




####
if len(sys.argv) != int(sys.argv[1])*6+5:
    sys.exit("Usage: python3 hemepureNMeshpipeline2.py NumMESHES MapORDER(as list of pairs 'outIDX,inIDX,out...'; first mesh to follow here is regarded as 0; set as '0,0' if doing single mesh only) {STLFNAME_A STLUNITS_A(e.g 1e-3 for mm) INLETPOSITIONS_A(X1,Y1,Z1;X2,Y2,Z2;..., (in quotes)) NUMINLETS_A NUMOUTLETS_A DXreq_A } NUMRANKS tauDesired; repeat part in curly braces for each mesh being analysed")

NUMRANKS = int(sys.argv[-2])
tauDes = float(sys.argv[-1])

numMeshes = int(sys.argv[1])
outletMeshes = [int(i) for i in (sys.argv[2]).split(',')][0::2]
inletMeshes = [int(i) for i in (sys.argv[2]).split(',')][1::2]

if numMeshes-1 != max(max(outletMeshes), max(inletMeshes)):
    sys.exit("Potential error in input -> number of meshes not match number of mapped meshes")

if len(outletMeshes) != len(inletMeshes):
    sys.exit("Potential error in mapping input -> number of inlets/outlets do not match")

inletDict={}
outletDict={}
iRadDict={}
oRadDict={}

for mesh in range(numMeshes):
    STLFNAME = sys.argv[3+mesh*6]
    STLUNITS = float(sys.argv[4+mesh*6])
    INLETS = [np.float_(iolet.split(",")) for iolet in (sys.argv[5+mesh*6]).split(";")]
    NUMINLETS = int(sys.argv[6+mesh*6])
    NUMOUTLETS = int(sys.argv[7+mesh*6])
    DXreq = np.float(sys.argv[8+mesh*6])

    ROOTNAME = os.path.splitext(os.path.basename(STLFNAME))[0]

    print("Values for mesh ", mesh)
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
    write_voxelizer_xml(xmlfname, DXreq, STLFNAME, inletpos0, outletpos0)

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
    
    if mesh == 0:
        shiftMaster = dx*shifts


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
    write_voxelizer_xml(xmlfname, DXreq, STLFNAME, inletposlist, outletposlist)

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
    write_heme_xml(mesh, tauDes, hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt, dx*shifts - shiftMaster)

    # Convert the voxelizer output into a hemelb gmy file
    execute("bash " + MAKEGMYMPIPATH + " fluidAndLinks"+str(mesh)+".dat " + gmyfname  + " " + str(NUMRANKS) + " " + str(VX2GMY_CHUNKSIZE) + "\n")
##################

# write mapping of outlets to inlets
print("Calculating and writing outlet/inlet mappings")

for i in range(len(outletMeshes)):
    write_dualMap(outletMeshes[i],inletMeshes[i])
    write_scalingMapFile(outletMeshes[i],inletMeshes[i])


## Create the velocity weights file - WARNING: CURRENTLY ASSUMES ONLY 1 INLET (not easy to fix...)
#inletsfname = ROOTNAME + ".inlets"
#execute(GMY2INLETSPATH + " " + gmyfname + " " + inletsfname + "\n")
#execute("python " + INFLOWPROFILEBUILDERPATH + " " + inletsfname + "\n")
#execute("cp out.weights.txt INLET0_VELOCITY.txt.weights.txt\n")
