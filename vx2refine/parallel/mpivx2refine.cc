/**
 *  Some technical notes about the HemeLB .gmy format:
 *
 * "There is a pleasure sure
 *  In being mad which none but
 *  Madmen know."
 *   -- John Dryden, The Spanish Friar, 1681
 *
 */

#include <stdio.h>
#include <cstring>
#include <string>
#include <stdint.h>
#include <math.h>
#include <rpc/rpc.h>
#include <zlib.h>
#include <unordered_map>
#include <map>
#include <vector>
#include <iostream>
#include <limits.h>
#include <vector>
#include <assert.h>
#include <mpi.h>

#define WARN_FOR_MEM_OVER_GB 1

#define DEBUG 0

#if DEBUG==1
#define DEBUGMSG(...) printf(__VA_ARGS__)
#else
#define DEBUGMSG(...)
#endif

using namespace std;


uint64_t sx[8] = {0,0,0,0,1,1,1,1};
uint64_t sy[8] = {0,0,1,1,0,0,1,1};
uint64_t sz[8] = {0,1,0,1,0,1,0,1};

int64_t linkMap[8][26]
{
	{ 0, 1, 1, 3, 4, 4, 3, 4, 4, 9,10,10,12,-1,12,-1,-1, 9,10,10,12,-1,-1,12,-1,-1},
	{ 1, 1, 2, 4, 4, 5, 4, 4, 5,10,10,11,-1,13,-1,-1,13,10,10,11,-1,-1,13,-1,-1,13},
	{ 3, 4, 4, 3, 4, 4, 6, 7, 7,12,-1,-1,12,-1,14,15,15,12,-1,-1,12,-1,-1,14,15,15},
	{ 4, 4, 5, 4, 4, 5, 7, 7, 8,-1,-1,13,-1,13,15,15,16,-1,-1,13,-1,-1,13,15,15,16},
	{ 9,10,10,12,-1,-1,12,-1,-1, 9,10,10,12,-1,12,-1,-1,17,18,18,20,21,21,20,21,21},
	{10,10,11,-1,-1,13,-1,-1,13,10,10,11,-1,13,-1,-1,13,18,18,19,21,21,22,21,21,22},
	{12,-1,-1,12,-1,-1,14,15,15,12,-1,-1,12,-1,14,15,15,20,21,21,20,21,21,23,24,24},
	{-1,-1,13,-1,-1,13,15,15,16,-1,-1,13,-1,13,15,15,16,21,21,22,21,21,22,24,24,25}
};

static const int64_t num_lints_per_record = 3 + 26 * 6;
static const int64_t size_of_record = num_lints_per_record * sizeof(uint64_t);
// convert long int to double
double plint_to_double(uint64_t i) {
    double result = 0;
    memcpy(&result, &i, sizeof(i));
    return result;
}

// convert double to long int
uint64_t double_to_plint(double i) {
    uint64_t result = 0;
    memcpy(&result, &i, sizeof(i));
    return result;
}
static const uint64_t ul0=0;
void coarsecell_to_subcells(uint64_t *data,uint64_t subcell, uint64_t *subcelldata) {
  
    subcelldata[0] = (uint64_t)( 2*data[0]+sx[subcell]); //x position
    subcelldata[1] = (uint64_t)( 2*data[1]+sy[subcell]); //y position
    subcelldata[2] = (uint64_t)( 2*data[2]+sz[subcell]); //z position
	
    for (uint64_t j=0; j<26; j++){
	    if (linkMap[subcell][j]>=0){
   	        subcelldata[6*j+3] = (uint64_t) data[6*linkMap[subcell][j]+3]; // linkType mapping
   	        
        	if (subcelldata[6*j+3]>0) { //is a wall
	            subcelldata[6*j+4] = (uint64_t) data[6*linkMap[subcell][j]+4]; // ConfigID mapping    
		    //Carry through coarse cell data, assuming that wall equivalently refines to grid
		    //this does move wall position slightly but remains outside original fluid cells
		    //i.e. wall asymptotically approaches edge of original fluid cell, max movement =1dx_orig
		    subcelldata[6*j+5] = (uint64_t) data[6*linkMap[subcell][j]+5]; //WallDistance
		    subcelldata[6*j+6] = (uint64_t) data[6*linkMap[subcell][j]+6]; //Nx
		    subcelldata[6*j+7] = (uint64_t) data[6*linkMap[subcell][j]+7]; //Ny
		    subcelldata[6*j+8] = (uint64_t) data[6*linkMap[subcell][j]+8]; //Nz
		} else {
                    for(uint64_t i=4;i<9;i++) subcelldata[6*j+i] = ul0; //nullify entries
                }
	    }
    }
}

int vx2refine(char *fname_vx, char *fnameRef_vx, uint64_t bsize, MPI_Comm comm)
{

	MPI_Offset total_file_size;
	MPI_Offset offset;
	int32_t error;
	int32_t this_rank;
	int32_t num_ranks;

	const int64_t num_lints_per_record = 3 + 26 * 6;
	const int64_t size_of_record = num_lints_per_record * sizeof(uint64_t);

	// get the current rank and total number of ranks
	MPI_Comm_rank(comm, &this_rank);
	MPI_Comm_size(comm, &num_ranks);

	// open vx output file for reading
	DEBUGMSG("Rank %d, opening %s...\n", this_rank, fname_vx);
	MPI_File vx_in;
	error = MPI_File_open(comm, fname_vx, MPI_MODE_RDONLY, MPI_INFO_NULL, &vx_in);
	if(error != MPI_SUCCESS) {
		fprintf(stderr, "ERROR: Could not open '%s' for reading. Error %d on rank %d.\n", fname_vx, error, this_rank);
		MPI_Finalize();
		exit(1);
	}

	// first check that there are an equal number of records in this file
	MPI_File_get_size(vx_in, &total_file_size);
	if(total_file_size % size_of_record != 0) {
		fprintf(stderr, "ERROR: File (size %lld) does not contain a whole number of records (of size %ld).\n", total_file_size, size_of_record);
		MPI_Finalize();
		exit(1);
	}

	bool store_records_on_rank0 = true;

	// work out what part of the file this rank should read
	uint64_t total_num_records = total_file_size/size_of_record;
	uint64_t num_records_per_rank = 0, num_records_on_rank0 = 0, num_records_on_rank1 = 0;

	// Decide how many records from the file will need to be read in by ranks 0, 1 and greater. If store_records_on_rank0 is false
	// (usually to stop rank0 going OOM on large block chunk sizes) then rank 1 will take the "excess" records.
	if(store_records_on_rank0) {
		num_records_per_rank = (uint64_t)round(total_num_records/(float)num_ranks + .5); // round up
		num_records_on_rank1 = num_records_per_rank;
		num_records_on_rank0 = total_num_records % num_records_per_rank;
	} else {
		num_records_per_rank = (uint64_t)round(total_num_records/(float)(num_ranks - 1) + .5); // round up
		num_records_on_rank1 = total_num_records % num_records_per_rank;
		num_records_on_rank0 = 0;
	}

	if(this_rank == 0)
		printf(
				"         Total number of records = %lu\n"
				"     Number of records on rank 0 = %lu\n"
				"     Number of records on rank 1 = %lu\n"
				"Number of records on other ranks = %lu\n", total_num_records, num_records_on_rank0, num_records_on_rank1, num_records_per_rank);

	// decide the offset at which to start reading the file for this rank
	if(this_rank == 0) {
		offset = 0;
	} else if(this_rank == 1) {
		offset = num_records_on_rank0 * size_of_record;
	} else {
		offset = (num_records_on_rank0 + num_records_on_rank1 + (this_rank - 2) * num_records_per_rank) * size_of_record;
	}

	// watch out for first rank (potentially) having a different number of records
	uint64_t num_records_on_this_rank;
	if(this_rank == 0) {
		num_records_on_this_rank = num_records_on_rank0;
	} else if(this_rank == 1) {
		num_records_on_this_rank = num_records_on_rank1;
	} else {
		num_records_on_this_rank = num_records_per_rank;
	}


    FILE *vxfile = fopen(fname_vx, "rb");
    if(vxfile == NULL) {
        fprintf(stderr, "Could not open file for reading.\n");
        return 1;
    }

    char fnameRefvx[32];
    snprintf(fnameRefvx, sizeof(char)*32, "refinedSubFile_%d.dat", this_rank);

    FILE *vxfileRef = fopen(fnameRefvx,"wb");
    // read bsize records
    uint64_t nread;
    uint64_t *record  = ( uint64_t *)calloc(  num_lints_per_record*bsize,sizeof(uint64_t));
    uint64_t *refined = ( uint64_t *)calloc(8*num_lints_per_record*bsize,sizeof(uint64_t));
    fprintf(stderr, "Blocksizes: %lu %lu.\n",size_of_record*bsize*8,size_of_record*bsize*8*8);

    uint64_t readRecords = 0;

    fseek(vxfile,offset,SEEK_SET);
	
    while( readRecords < (num_records_on_this_rank - bsize)){
	nread=fread(record, 1, size_of_record*bsize , vxfile);
	readRecords += bsize;
        //memset (refined,8*num_lints_per_record*bsize*sizeof(uint64_t),0);
	for (uint64_t i=0; i<bsize; i++) {
	   for (uint64_t j=0; j<8; j++) {
	       coarsecell_to_subcells(record+num_lints_per_record*i,j,refined+num_lints_per_record*(i*8+j));
	   }
        }
        fwrite(refined,1, 8* size_of_record *bsize, vxfileRef);
	fflush(vxfileRef);
    }
    
    //get the last bit
    fread(record, 1, size_of_record*bsize , vxfile);
    // remainder loop
    nread=num_records_on_this_rank - readRecords;
    //memset (refined,8*num_lints_per_record*nread*sizeof(uint64_t),0);
    for (uint64_t i=0; i<nread; i++) {
	for (uint64_t j=0; j<8; j++) {
	    coarsecell_to_subcells(record+num_lints_per_record*i,j,refined+num_lints_per_record*(i*8+j));
	}
    }
    fwrite(refined,1, 8* size_of_record*nread, vxfileRef);
    fflush(vxfileRef);

    fclose(vxfile);
    fclose(vxfileRef);
    free(record);
    free(refined);

    return 0;
}

int main(int argc, char **argv)
{
	assert (sizeof(uint64_t) == sizeof(double));

	if (argc != 3 && argc != 4) {
		fprintf(stderr, "Usage: vx2refine file.dat fileREF1.dat [BLOCKCHUNKSIZE].\n");
		return 1;
	}

	char *fname_vx = argv[1];
	char *fnameRef_vx = argv[2];
	uint64_t bsize = 1;
	if(argc == 4) {
		bsize = atoi(argv[3]);
	}

	// set up MPI
	MPI_Init(NULL, NULL);

	uint32_t ret = vx2refine(fname_vx, fnameRef_vx, bsize, MPI_COMM_WORLD);
	if (ret != 0) {
		fprintf(stderr, "Error ocurred.\n");
		return 1;
	}

	MPI_Finalize();

	return 0;
}
