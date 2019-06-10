/**
 *
 *  Some technical notes about the HemeLB .gmy format:
 *
 * "There is a pleasure sure
 *  In being mad which none but
 *  Madmen know."
 *   -- John Dryden, The Spanish Friar, 1681
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <rpc/rpc.h> 
#include <zlib.h>
#include <iostream>

#define DEBUG 0 

#if DEBUG==1
#define DEBUGMSG(...) printf(__VA_ARGS__)
#else
#define DEBUGMSG(...)
#endif

using namespace std;

int64_t xyz_to_i(uint32_t x, uint32_t y, uint32_t z, uint32_t maxX, uint32_t maxY,  uint32_t maxZ)
{
  if(x < 0 || x >= maxX || y < 0 || y >= maxY || z < 0 || z >=maxZ) {
    return -1;
  }
  return int64_t(x)*int64_t(maxY)*int64_t(maxZ) + int64_t(y)*int64_t(maxZ) + int64_t(z);
}

void unique_to_xyz(int64_t unique, uint32_t *X, uint32_t *Y, uint32_t *Z, uint32_t maxX, uint32_t maxY, uint32_t maxZ)
{
  *Z = unique%maxZ;
  int64_t temp = unique/maxZ;
  *Y = temp%maxY;
  *X = temp/maxY;
}

int write_uint(FILE* file, uint32_t value) {
  int rv = fwrite(&value,4,1,file);

  return rv;
}
int read_uint(FILE* file, uint32_t* value) {
  int rv = fread(value,4,1,file);

  return rv;
}


int gmy2inlets(char *fname_gmy, char *fname_inlets)
{
  // gmy file handle
  FILE *gmyFile;
  
  // xdr stdio stream
  XDR xdrs;
  // xdr block data / memory stream
  XDR xdrbs;

  // iterator
  uint32_t i;
  // return value
  uint32_t ret;

  // Declarations for the header information
  uint32_t hlbMagic;
  uint32_t gmyMagic;
  uint32_t gmyVersion;
  uint32_t xBlockCount;
  uint32_t yBlockCount;
  uint32_t zBlockCount;
  uint32_t blockEdgeLen;
  uint32_t padding;

  // Data calculated from header information
  uint32_t totalBlockCount;
  uint32_t maxSitesPerBlock;

  // Declarations for the block information header
  uint32_t *blockSiteCount;
  uint32_t *blockCompressedLen;
  uint32_t *blockUncompressedLen;
  // uint32_t *blockDataOffset;

  // Declarations for zlib decompression
  z_stream strm;
  unsigned char *compressedBuffer; 
  unsigned char *decompressedBuffer; 

  // Declarations for site information 
  uint32_t siteIsSimulated;
  uint32_t linkType;
  uint32_t siteType;
  float wallDistance;
  uint32_t configEntry;
  uint32_t hasWallNormal;
  float xWallNormal, yWallNormal, zWallNormal;

  uint32_t blockX = 0, blockY = 0, blockZ = 0;
  uint32_t localX = 0, localY = 0, localZ = 0;
 // uint32_t globalX = 0, globalY = 0, globalZ = 0;

  // Stuff needed for hgb
  uint32_t globalSiteCounter = 0;

  // Allocations to XDR_DECODE from file stream
  gmyFile = fopen(fname_gmy, "r");
  xdrstdio_create(&xdrs, gmyFile, XDR_DECODE);
  i = 0;

  // Open output inlets file
  FILE *outfile = NULL;
  if((outfile = fopen(fname_inlets, "w")) == NULL) {
	fprintf(stderr, "Could not open %s for writing. Aborting.\n", fname_inlets);
	return 100;
  }

  // Allocations for simple zlib decompression
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;

  // Reading gmy header information
  xdr_u_int(&xdrs, &hlbMagic);
  xdr_u_int(&xdrs, &gmyMagic);
  xdr_u_int(&xdrs, &gmyVersion);
  xdr_u_int(&xdrs, &xBlockCount);
  xdr_u_int(&xdrs, &yBlockCount);
  xdr_u_int(&xdrs, &zBlockCount);
  xdr_u_int(&xdrs, &blockEdgeLen);
  xdr_u_int(&xdrs, &padding);

  DEBUGMSG("hlbMagic = %u\n", hlbMagic);
  DEBUGMSG("gmyMagic = %u\n", gmyMagic);
  DEBUGMSG("gmyVersion = %u\n", gmyVersion);
  DEBUGMSG("xBlockCount = %u\n", xBlockCount);
  DEBUGMSG("yBlockCount = %u\n", yBlockCount);
  DEBUGMSG("zBlockCount = %u\n", zBlockCount);
  DEBUGMSG("blockEdgeLen = %u\n", blockEdgeLen);
  DEBUGMSG("padding = %u\n", padding);
  
  // Calculating the number of blocks filling the complete volume
  totalBlockCount = xBlockCount * yBlockCount * zBlockCount;
  // Calculating the maximum number of sites per block
  maxSitesPerBlock = blockEdgeLen * blockEdgeLen * blockEdgeLen;
  
  // Allocating arrays to read block header information
  blockSiteCount = new uint32_t[totalBlockCount];
  blockCompressedLen = new uint32_t[totalBlockCount];
  blockUncompressedLen = new uint32_t[totalBlockCount];
  
  // blockDataOffset = malloc(totalBlockCount * 3 * sizeof(uint32_t));  
  // blockDataOffset[0] = (8 + totalBlockCount * 3 ) * sizeof(uint32_t);

  // Reading block header information to arrays
  uint32_t num_inlets_found = 0;
  uint32_t num_entries = 0;
  while (i<totalBlockCount) {

    // Number of sites actually stored for the block
    xdr_u_int(&xdrs, &blockSiteCount[i]);
    // Compressed size of block data in bytes
    xdr_u_int(&xdrs, &blockCompressedLen[i]);
    // Uncompressed size of block data in bytes
    xdr_u_int(&xdrs, &blockUncompressedLen[i]);

    DEBUGMSG("Block %u:\n", i);
    DEBUGMSG("blockSiteCount = %u\n", blockSiteCount[i]);
    DEBUGMSG("blockCompressedLen = %u\n", blockCompressedLen[i]);
    DEBUGMSG("blockUncompressedLen = %u\n", blockUncompressedLen[i]);
 
    num_entries += blockSiteCount[i];

    // printf("i: %d %d %d %d\n",i,blockSiteCount[i],blockCompressedLen[i],blockUncompressedLen[i]);
    // Block data offset can be used to skip directly to blocks / parallel access
    // blockDataOffset[i+1] = blockCompressedLen[i] + blockDataOffset[i];
    i++;
  }

  // Reiterate blocks, this time read actual data sets
  for (i=0;i<totalBlockCount;i++) {
    DEBUGMSG("Block %u:\n", i);
    // There is only block data available if there is at least one lattice site on the block
    if (blockSiteCount[i] == 0) {
      DEBUGMSG("No lattice sites on this block.\n");
    } else {
      // Allocate respective buffers
      compressedBuffer =  new unsigned char[blockCompressedLen[i]];
      decompressedBuffer = new unsigned char[blockUncompressedLen[i]];
      // Read compressed block from the file
      ret = fread(compressedBuffer, blockCompressedLen[i], 1, gmyFile);
      // Initialise zlib and check if it is OK
      ret = inflateInit(&strm);
      if (ret != Z_OK)
        printf("Decompression init error for block\n"); 

      // Specify number of bytes to decompress (inflate in zlib-speech)
      strm.avail_in = blockCompressedLen[i];
      // Specify compressed input buffer
      strm.next_in = compressedBuffer;

      // Specify number of uncompressed bytes expected after decompression (inflate in zlib-speech)
      // (This is normally used to to read chunks to/from a stream)
      strm.avail_out = blockUncompressedLen[i];
      // Specify compressed output buffer
      strm.next_out = decompressedBuffer;

      // Run zlib decompression and check if it is OK
      ret = inflate(&strm, Z_FINISH);
      if (ret != Z_STREAM_END)
        printf("Decompression error for block\n");

      // Finalise/end zlib decompression and check if it is OK
      ret = inflateEnd(&strm);
      if (ret != Z_OK)
        printf("Decompression end error for block\n");

      delete(compressedBuffer);
          
      // Create a new XDR_DECODE stream on decompressed buffer
      xdrmem_create(&xdrbs, (char *) decompressedBuffer, blockUncompressedLen[i],XDR_DECODE);

      // Loop through all sites on block
      uint32_t site = 0;
      //      while(site < blockSiteCount[i]) {
      //	for (site = 0;site < blockSiteCount[i];site++) {
      localX = 0;
      localY = 0;
      localZ = 0;

      for (site = 0;site < maxSitesPerBlock; site++) {
        DEBUGMSG("local site index (on block): %u\n", site);
        // Check if inside or outside (aka solid) of the simulation domain
        xdr_u_int(&xdrbs, &siteIsSimulated);
        DEBUGMSG("siteIsSimulated=%u\n", siteIsSimulated);
        if (siteIsSimulated == 1) {
          DEBUGMSG("globalSiteCounter: %u\n", globalSiteCounter);
	if (globalSiteCounter%100000 == 0) {
	  printf("Read site %u of %u\n",globalSiteCounter,num_entries);

	}
          // Default type is FLUID = 0
          siteType = 0;

          uint32_t link = 0;
          for (link=0;link<26;link++) { // Only consider 26 links (because there is no link to itself)
            // Check type for every link to our neighbours
            xdr_u_int(&xdrbs, &linkType);
            DEBUGMSG("Link %u: linkType = %u\n", link, linkType);
            switch (linkType) {
            case 0 :
              // linkType = FLUID
              DEBUGMSG("Link %u: FLUID - no further data.\n", link);
              break;
            case 1 : 
              // linkType = WALL
              // Read float of distance from nearest obstacle
              xdr_float(&xdrbs, &wallDistance);
              DEBUGMSG("Link %u: WALL - wallDistance=%f\n", link, wallDistance);
              switch (siteType) {
                // Check for combinations of linkTypes to specify combined siteTypes
              case 0 : // was FLUID
                siteType = 1; // is now WALL
                break;
              case 1 : // was WALL
                siteType = 1; // stays WALL
                break;
              case 2 : // was INLET
                siteType = 4; // is now INLET+WALL
                break;
              case 3 : // was OUTLET
                siteType = 5; // is now OUTLET+WALL
                break;
              case 4 : // was INLET+WALL
                siteType = 4; // stays INLET+WALL
                break;
              case 5 : // was OUTLET+WALL
                siteType = 5; // stays OUTLET+WALL
                break;
              }
              break;
            case 2 :
              // linkType = INLET
              // Read float of distance from nearest obstacle
              xdr_u_int(&xdrbs, &configEntry);
              // Read float of distance from nearest obstacle
              xdr_float(&xdrbs, &wallDistance);
              DEBUGMSG("Link %u: INLET - configEntry=%u\n", link, configEntry);
              DEBUGMSG("Link %u: INLET - wallDistance=%f\n", link, wallDistance);
              switch (siteType) {
              case 0 : // was FLUID
                siteType = 2; // is now INLET
                break;
              case 1 : // was WALL
                siteType = 4;  // is now INLET+WALL
                break;
              case 2 : // was INLET
                siteType = 2; // stays INLET
                break;
              case 3 : // was OUTLET
                siteType = -2; // INLET+OUTLET = Kaputt
                break;
              case 4 : // was INLET+WALL
                siteType = 4; // stays INLET+WALL
                break;
              case 5 : // was OUTLET+WALL
                siteType = -2; // INLET+OUTLET = Kaputt
                break;
              }
              break;
            case 3 :
              // linkType = OUTLET
              // Read index of associated configuration in config.xml
              xdr_u_int(&xdrbs, &configEntry);
             // Read float of distance from nearest obstacle
              xdr_float(&xdrbs, &wallDistance);
              DEBUGMSG("Link %u: OUTLET - configEntry=%u\n", link, configEntry);
              DEBUGMSG("Link %u: OUTLET - wallDistance=%f\n", link, wallDistance);
              switch (siteType) {
              case 0 : // was FLUID
                siteType = 3; // is now OUTLET
                break;
              case 1 : // was WALL
                siteType = 5; // is now OUTLET+WALL
                break;
              case 2 : // was INLET
                siteType = -2; // INLET+OUTLET = Kaputt
                break;
              case 3 : // was OUTLET
                siteType = 3; // stays OUTLET
                break;
              case 4 : // was INLET+WALL
                siteType = -2; // INLET+OUTLET = Kaputt
                break;
              case 5 : // was OUTLET+WALL
                siteType = 5; // stays OUTLET+WALL
                break;              
              }
              break;

            default:
              printf("Unrecognised linkType %u\n", linkType);
              exit(1);
            }
          } // link loop

	  uint32_t globalX = 0, globalY = 0, globalZ = 0;
          globalX = blockX * blockEdgeLen + localX;
          globalY = blockY * blockEdgeLen + localY;          
          globalZ = blockZ * blockEdgeLen + localZ;

          // Check if there are wall normal coordinates to be read
          xdr_u_int(&xdrbs, &hasWallNormal);
          DEBUGMSG("hasWallNormal=%u\n", hasWallNormal);
          if (hasWallNormal == 1) {
            // Read wall normal coordinates as separate floats
            xdr_float(&xdrbs, &xWallNormal);
            xdr_float(&xdrbs, &yWallNormal);
            xdr_float(&xdrbs, &zWallNormal);
            DEBUGMSG("xWallNormal=%f\n", xWallNormal);
            DEBUGMSG("yWallNormal=%f\n", yWallNormal);
            DEBUGMSG("zWallNormal=%f\n", zWallNormal);
          } else if (hasWallNormal != 0) {
            printf("hasWallNormal has some messed up value, bro. (%u)\n", hasWallNormal);
            exit(1);
          }


	  // Check if site is inlet
	  if(siteType == 2 || siteType == 4) {
		  num_inlets_found += 1;
	          fprintf(outfile, "%d %d %d %d\n", globalX, globalY, globalZ, siteType);
	  }

          globalSiteCounter += 1;
          //          site += 1;
        } else {
          // a solid site does not need a site type
          siteType = -1;
        }
       
        localZ += 1;
        if(localZ == blockEdgeLen) {
          localZ = 0;
          localY += 1;
          if(localY == blockEdgeLen) {
            localY = 0;
            localX += 1;
          }
        }
        
      } // site loop      
      // Destroy the XDR_DECODE stream on decompressedBuffer
      xdr_destroy(&xdrbs);     
      delete(decompressedBuffer); 
    } // if block not empty

    blockZ += 1;
    if(blockZ == zBlockCount) {
      blockZ = 0;
      blockY += 1;
      if(blockY == yBlockCount) {
        blockY = 0;
        blockX += 1;
      }
    }
    
  }
    
  printf("Found %d inlet sites.\n", num_inlets_found);

  // Free arrays used to store block header information
  delete(blockSiteCount);
  delete(blockCompressedLen);
  delete(blockUncompressedLen);

  // free(blockDataOffset);

  // Destroy the XDR_DECODE stream on the gmy file
  xdr_destroy(&xdrs);

  // Close the gmy file handle
  fclose(gmyFile);
  fclose(outfile);

  return 1;
}



int main(int argc, char **argv)
{
  if (argc == 3) {

    char *fname_gmy = argv[1];
    char *fname_inlets = argv[2];

    int ret = gmy2inlets(fname_gmy, fname_inlets);
    if (ret != 1) {
      return 1;
    }
    return 0;
  }
  /* otherwise, report usage */
  else {
    fputs("Usage: gmy2inlets infile.gmy outfile.inlets\n", stderr);
    return 1;
  }
}

