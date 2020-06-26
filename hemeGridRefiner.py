import os, sys
import numpy as np

VOXELIZERPATH = "/cs/heme/HemePure_JM/HemePure_tools/voxelizer/source/voxelizer3"#_file_ARCHER"
MAKEGMYMPIPATH = "/cs/heme/HemePure_JM/HemePure_tools/vx2gmy/make_gmy_MPI.sh"
GMY2INLETSPATH = "/cs/heme/HemePure_JM/HemePure_tools/gmy2inlets/gmy2inlets"
VXREFINERPATH = "/cs/heme/HemePure_JM/HemePure_tools/vx2refine/vx2refine"
VXPARAVIEWPATH = "/cs/heme/HemePure_JM/HemePure_tools/vx2paraview/vx2paraview"
#INFLOWPROFILEBUILDERPATH = "../inflow-profile-builder/inflow.py"

VX2GMY_CHUNKSIZE = 2000

def execute(command):
        print("Executing: " + command)
        r = os.system(command)
        if r != 0:
                sys.exit("Command failed.")

def write_heme_xml(tauValue, hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt):
    xml =  "<?xml version=\"1.0\"?>\n"
    xml += "<hemelbsettings version=\"3\">\n"
    xml += "  <simulation>\n"
    xml += "    <step_length units=\"s\" value= " + str((tauValue - 0.5)*gmy_resolution*gmy_resolution/(12e-6)) + "/>\n"
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


if len(sys.argv) != 6:
    sys.exit("Usage: python3 hemeGridRefine.py STLFNAME NumRefineCycles DXorig NumRanks tauDesired")


tauDes = float(sys.argv[-1])
NUMRANKS = np.int(sys.argv[-2])

STLFNAME = sys.argv[1]
ROOTNAME = os.path.splitext(os.path.basename(STLFNAME))[0]
REFINES = np.int(sys.argv[2])
DXorig = np.float(sys.argv[3])

execute(VXREFINERPATH + " fluidAndLinks.dat" + " fluidAndLinks_Ref1.dat") 
for cycles in range(1,REFINES):
    execute(VXREFINERPATH + " fluidAndLinks_Ref" + str(cycles) + ".dat" + " fluidAndLinks_Ref" + str(cycles+1) + ".dat")
    execute("rm  fluidAndLinks_Ref" + str(cycles) + ".dat")

    with open("iolets_block_inputxml.txt", "r") as f:
        lines = f.readlines()
        j=0
        for i in lines:
            i = i.split()

            if i[0] == '<normal':
                normal = np.float_(i[2][8:-4].split(","))
            if i[0] == '<position':
                centre = np.float_(i[2][8:-4].split(","))
                centre = np.around(2*centre + 0.5*np.ones(np.shape(centre)) - normal, decimals=2)
                lines[j] = '      <position units="lattice" value="(' \
                           + str(centre[0]) + ',' + str(centre[1]) + ',' + str(centre[2]) + ')"/>\n'
            j = j+1

    with open("iolets_block_inputxml.txt", 'w') as f:
        f.writelines( lines )

# Get the inlets and outlets xml blocks (for the hemelb input xml) output by the voxelizer
with open("iolets_block_inputxml.txt", "r") as ioletsblockfile:
    ioletsblocktxt = ioletsblockfile.read()

# Write the hemelb input.xml file
hemexmlfname = "inputREF.xml"
gmyfname = ROOTNAME + ".gmy"
gmy_resolution = DXorig/2**REFINES 

write_heme_xml(tauDes, hemexmlfname, gmyfname, gmy_resolution, ioletsblocktxt) 

# Convert the refined voxelizer output into a hemelb gmy file
execute("bash " + MAKEGMYMPIPATH + " fluidAndLinks_Ref" + str(REFINES) + ".dat " + gmyfname  + " " + str(NUMRANKS) + " " + str(VX2GMY_CHUNKSIZE) + "\n")

# Convert refined voxelizer output into a paraview readable form to allow geometry check
execute(VXPARAVIEWPATH + " fluidAndLinks_Ref" + str(REFINES) + ".dat > " + ROOTNAME + "RefCHECK.txt") 

## Create the velocity weights file 
#inletsfname = ROOTNAME + ".inlets"
#execute(GMY2INLETSPATH + " " + gmyfname + " " + inletsfname + "\n")
#execute("python3 " + INFLOWPROFILEBUILDERPATH + " " + inletsfname + "\n")
#
#for ilet in range(0,NUMINLETS): 
#    execute("cp out" + str(ilet) + ".weights.txt INLET" + str(ilet) + "_VELOCITY.txt.weights.txt\n") 

