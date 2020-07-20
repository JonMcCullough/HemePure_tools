#include <cstring>
#include <string>
#include <stdlib.h> 
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
#include <cstdio> 

#define DEBUG 0

#if DEBUG==1
#define DEBUGMSG(...) printf(__VA_ARGS__)
#else
#define DEBUGMSG(...)
#endif

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

int main(int argc, char **argv)
{
    assert (sizeof(uint64_t) == sizeof(double));

    if (argc < 3 ||  argc>4) {
	fprintf(stderr, "Usage: vx2ascii file.vx file_refined.vx [blocksize]\n");
	return 1;
    }

    char *fname_vx = argv[1];
    char *fnameRef_vx = argv[2];
    uint64_t bsize=1;
    if (argc==4){
      bsize = atoi(argv[3]);
    }

    FILE *vxfile = fopen(fname_vx, "rb");
    if(vxfile == NULL) {
        fprintf(stderr, "Could not open file for reading.\n");
        return 1;
    }

    FILE *vxfileRef = fopen(fnameRef_vx,"wb");
    // read bsize records
    uint64_t nread;
    uint64_t *record  = ( uint64_t *)calloc(  num_lints_per_record*bsize,sizeof(uint64_t));
    uint64_t *refined = ( uint64_t *)calloc(8*num_lints_per_record*bsize,sizeof(uint64_t));
    fprintf(stderr, "Blocksizes: %lu %lu.\n",size_of_record*bsize*8,size_of_record*bsize*8*8);
    while( (nread=fread(record, 1, size_of_record*bsize , vxfile)) == size_of_record*bsize ) {
        //memset (refined,8*num_lints_per_record*bsize*sizeof(uint64_t),0);
	for (uint64_t i=0; i<bsize; i++) {
	   for (uint64_t j=0; j<8; j++) {
	       coarsecell_to_subcells(record+num_lints_per_record*i,j,refined+num_lints_per_record*(i*8+j));
	   }
        }
        fwrite(refined,1, 8* size_of_record *bsize, vxfileRef);
	fflush(vxfileRef);
    }
    // remainder loop
    nread/=size_of_record;
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
