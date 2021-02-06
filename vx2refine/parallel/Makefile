MPICC=mpiCC#/opt/cray/craype/2.5.8/bin/CC

mpivx2refine: mpivx2refine.cc
	$(MPICC) -O3 -Wall -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -std=c++11 mpivx2refine.cc -lz -o mpivx2refine
debug: mpivx2refine.cc
	$(MPICC) -O0 -g -Wall -D_LARGEFILE_SOURCE=1 -D_FILE_OFFSET_BITS=64 -std=c++11 mpivx2refine.cc -lz -o mpivx2refine
