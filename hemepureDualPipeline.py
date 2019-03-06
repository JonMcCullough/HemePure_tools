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

## JM stuff...
def write_dualMap(Aout, Vin):
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
    
    with open("outlets_radius_A.txt","r") as oletsA:
        oLETSA = [line.rstrip('\n') for line in oletsA]
        oLETSA = [float(i.split(',')[1]) for i in oLETSA]
        oletsA.close()

    with open("inlets_radius_B.txt","r") as iletsB:
        iLETSB = [line.rstrip('\n') for line in iletsB]
        iLETSB = [float(i.split(',')[1]) for i in iLETSB]
        iletsB.close()

    mapp = "Map Indices: Aout (index, radius) then Paired Vin(s) (index, distance, radius) \n"
    for o in range(nA):
        mapp += "(" + str(o+1) + "," + str(oLETSA[o]) + ")"
        for i in range(nI):
            if VinPaired[i] == o: mapp += ", ( " + str(i+1) + "," + str(dMinPaired[i]) + "," + str(iLETSB[i]) + ")"
        mapp += "\n"
    
    with open("mapAtoB.txt", "w") as outMap:
        outMap.write(mapp)

def NetworkCalc(r0,r1,L,ratio,scale):
    # Computation of network of vessels between outlet and an inlet
    # r0 = radius of outlet
    # r1 = radius of inlet being connected
    # ratio = rate at which daughter branches are narrowed compared to parent
    # scale = length of each branch segment - L = scale*radius_atlevel

    if L>scale*(r0+r1)/(1-ratio):
        print("Length check failed")
        exit()

    Lrem = L - scale*(r0+r1)
    branches_a = [r0]
    branches_v = [r1]

    d=0
    while Lrem >= 0.0:
        d = d+1
        dA = ratio*branches_a[-1]
        dV = ratio*branches_v[-1]
        print( "daughter level", d)
        print( "dA = ", dA, "; dV = ", dV)
        Lrem = Lrem - scale*(dA+dV)
        print( Lrem)

        if Lrem < 0.0:
            Lrem = Lrem + scale*(dA+dV)
            d = d-1
            print( "Extra straight connection of ", Lrem, "at daughter layer", d)
            break

        if (max(dA,dV)<1e-6):
            print( "Effective Connection Reached")
            break
        
        branches_a.append(dA)
        branches_v.append(dV)

    print( branches_a)
    print( branches_v)

    rla = 0; la = 0;
    j=0
    for i in branches_a:
        rla = rla + scale*i*i*2**j
        la = la + scale*i*2**j
        j = j+1
        
    rla = rla+Lrem*ratio*(branches_a[-1] + branches_v[-1])*2**(j-1)
    la = la+ratio*Lrem*2**(j-1)
    print( "average artery radius = ", rla/la)

    rlv = 0; lv = 0; j=0
    for i in branches_v:
        rlv = rlv + scale*i*i*2**j
        lv = lv + scale*i*2**j
        j = j+1
        
    rlv = rlv+Lrem*ratio*(branches_a[-1] + branches_v[-1])*2**(j-1)
    lv = lv+ratio*Lrem*2**(j-1)
    print( "average vein radius = ", rlv/lv)

    effR = (rla + rlv)/(la + lv)
    effL = la + lv

    if d == 0:
        effR = 0.5*(branches_a[0] + branches_v[0])
        effL = L
    
    print( "average vessel radius = ", effR)
    print( "effective length = ", effL)

    print( "effective volume = ", np.pi*effR*effR*effL)
    print( "\n")
    return [effR, effL, np.pi*effR*effR*effL]
#

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
execute("mv inlets_radius.txt inlets_radius_A.txt")
execute("mv outlets_radius.txt outlets_radius_A.txt")

execute("mpirun -np " + str(NUMRANKS) + " " + VOXELIZERPATH + " " + xmlfname_B + "\n")
execute("cat fluidAndLinks_*.dat > fluidAndLinksB.dat && rm fluidAndLinks_*.plb && rm fluidAndLinks_*.dat")
execute("mv iolets_block_inputxml.txt iolets_block_inputxml_B.txt")
execute("mv inlets_radius.txt inlets_radius_B.txt")
execute("mv outlets_radius.txt outlets_radius_B.txt")

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

# compute boundary relationship values based on geometric values and outlet pressure (P0) and flow rate (Q0)
mu = 0.1666666 #3.5e-3 = Physical units, lattice units = 0.1666666(?)
frict = 1.3 #1.5 for turbulent, 1.0 for laminar
scale = 50 #typical value
ratio = 0.5 # Murray Law splitting

flow = "Flow Rate scale factors: \n"
pressure = "dP scale factors: \n"


with open("mapAtoB.txt") as theMap:
    next(theMap)
    for line in theMap:

        line = line.strip("\n").strip("(").strip(")").split("), (")
        outletData = [float(x) for x in line[0].split(",")]
        r0 = outletData[1]
        outletIDX = outletData[0]

        rInlets = np.zeros(len(line)-1)
        L_Inlets = np.zeros(len(line)-1)
        inletIDX = np.zeros(len(line)-1)
        for out in range(1,len(line)):
            inletData = [float(x) for x in line[out].split(",")]
            rInlets[out-1] = inletData[2]
            L_Inlets[out-1] = inletData[1]
            inletIDX[out-1] = inletData[0]

        A = np.zeros(len(rInlets))
        k = np.zeros(len(rInlets))
        v = np.zeros(len(rInlets))

        for i in range(len(rInlets)):
            print("Branch from outlet", outletIDX, "to inlet", inletIDX[i])
            r,e1,vol = NetworkCalc(r0,rInlets[i], L_Inlets[i],ratio,scale)

            A[i] = np.pi*r*r
            k[i] = (8*np.pi * mu)/(A[i]*A[i])
            v[i] = vol


        q = np.divide(v,np.sum(v))
        dP = frict*np.multiply(np.multiply(q,L_Inlets),k)

        print("WARNING - will need to run dP check in hemeLB to ensure ok")
        print( "Factors")
        print( "Q = ",q,"*Q0 ")
        print( "dP = ", dP,"*Q0 [L.U.]")
        
        for i in range(len(rInlets)):
            flow += str(int(outletIDX)) + "," + str(int(inletIDX[i])) + "," + str(q[i]) + "\n"
            pressure += str(int(outletIDX)) + "," + str(int(inletIDX[i])) + "," + str(dP[i]) + "\n"

with open("flowFactors.txt", "w") as factors:
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
