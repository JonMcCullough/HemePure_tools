# HemePure_tools
- Repository contains tools for HemeLB geometry generation and stuff.

- requires geometry in *.stl format with in/outlets open (think pipe) and inward oriented normals
- set-up to use Python3

1) clone repository into desired locations
2) modify paths at top of hemeLBpreproc.py (VOXELIZERPATH and MAKEGMYMPIPATH) and vx2gmy/make_gmy_MPI.sh (VX2GMYPATH)
3) run 'make' in voxelizer/source/, gmy2inlets/ and vx2gmy/. It appears voxelizer/source is sensitive to gcc version (has been successful with gcc/4.8 and gcc/5.5 but not gcc/7.4)
4) navigate to folder with *.stl file and run script with following instruction:
	python3 PATHTO/hemeLBpreprocSingle.py STLFNAME STLUNITS(e.g 1e-3 for mm) INLETPOSITIONS(X1,Y1,Z1;X2,Y2,Z2;..., (in quotes)) NUMINLETS NUMOUTLETS DXreq NUMRANKS tauDesired

	STLFNAME - name of *.stl file
	STLUNITS - units of *.stl file (e.g 1e-3 for mm)
	INLETPOSITIONS - approximate location of inlets in physical units formatted as: "X1,Y1,Z1;X2,Y2,Z2;...", hemeLBPreProcSingle.py expects this to be in a text file
	NUMINLETS - number of inlets to geometry
	NUMOUTLETS - number of outlets of geometry
	DXreq - required lattice spacing (in m)
	NUMRANKS - number of ranks to run parallel voxeliztion on
        tauDesired - desired value of LBM relaxation time, tau. Must be 0.5<tau and ideally between 0.9-1.0 if able. Governed by DXreq and HemeLB blood viscosity (0.004 Pa.s)

e.g. for pipe case in PipeTest/
	python3 ../hemeLBpreprocSingle.py pipe.stl 1.0 inlets.txt 1 1 0.1 2 1.0; where inlets.txt contains 0,0,0

5) once complete open the created input.xml file to check/edit parameters to desired values. In particular:
	<steps units="lattice" value="CHANGE"/> - modify CHANGE to desired simulation length
	<propertyoutput file="inlet.dat" period="100"> - modify period to desired output rate, also for outlet.dat and whole.dat (remove sections if desired)
	<mean units="mmHg" value="0.0"/> - modify to desired presssure value at inlet/outlet; other parameters here allow for an oscilliatory pressure profile. N.B. Pressure here refers to a pressure difference from the static LBM pressure (1/3 in lattice units, the difference needs to be small compared to this), see also initial condition

