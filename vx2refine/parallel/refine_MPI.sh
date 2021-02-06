# !/bin/bash

VX2GMYPATH=/net/storeptr1/heme/HemePure_JM/HemePure_tools/vx2refine/parallel

# Check that enough arguments have been given
if [[ $# -ne 4 ]]; then
        echo "Usage: bash make_gmy_MPI.sh INPUT.VOXEL OUTPUT.GMY NUMRANKS BLOCKCHUNKSIZE"
        exit 1
fi

mpirun -np $3 $VX2GMYPATH/mpivx2refine $1 __make_gmy_MPI_test.gmy $4
cat refinedSubFile_* > $2
rm refinedSubFile_*
