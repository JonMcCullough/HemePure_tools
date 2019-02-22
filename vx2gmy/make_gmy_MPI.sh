# !/bin/bash

VX2GMYPATH=/home/robinr/hemevox/HemePure/src_tools/vx2gmy

# Check that enough arguments have been given
if [[ $# -ne 4 ]]; then
        echo "Usage: bash make_gmy_MPI.sh INPUT.VOXEL OUTPUT.GMY NUMRANKS BLOCKCHUNKSIZE"
        exit 1
fi

mpirun -np $3 $VX2GMYPATH/mpivx2gmy $1 __make_gmy_MPI_test.gmy $4 && cat __make_gmy_MPI_test.gmy __make_gmy_MPI_test.gmy_blockdata_temp_ > $2
#aprun -n $3 -N 1 ./mpivx2gmy $1 __make_gmy_MPI_test.gmy $4 && cat __make_gmy_MPI_test.gmy __make_gmy_MPI_test.gmy_blockdata_temp_ > $2
