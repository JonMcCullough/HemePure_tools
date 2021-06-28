/**
 *
 *  Some technical notes about the HemeLB .gmy format:
 *
 * "There is a pleasure sure
 *  In being mad which none but madmen know."
 *   -- John Dryden, The Spanish Friar, 1681
 *
 *
 */

#include <iostream>
#include <math.h>
#include <rpc/rpc.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <tr1/unordered_map>
#include <zlib.h>

#define DEBUG 0

#if DEBUG==1
#define DEBUGMSG(...) printf(__VA_ARGS__)
#else
#define DEBUGMSG(...)
#endif

using namespace std;

int gmy2gmyplus(char *fname_gmy, char *fname_gmyplus, uint32_t weight_FLUID, uint32_t weight_WALL, uint32_t weight_INLET, uint32_t weight_OUTLET, uint32_t weight_INLETWALL, uint32_t weight_OUTLETWALL)
{
	// gmy file handle
	FILE *gmyFile = NULL;

	// xdr stdio stream
	XDR xdrs;
	// xdr block data / memory stream
	XDR xdrbs;
	// xdr stdio stream to write output file
	XDR xdrsw;

	// iterator
	uint32_t i;
	// return value
	uint32_t ret;

	// declarations for the header information
	uint32_t hlbMagic;
	uint32_t gmyMagic;
	uint32_t gmyVersion;
	uint32_t xBlockCount;
	uint32_t yBlockCount;
	uint32_t zBlockCount;
	uint32_t blockEdgeLen;
	uint32_t padding;

	// data calculated from header information
	uint32_t totalBlockCount;
	uint32_t maxSitesPerBlock;

	// declarations for the block information header
	uint32_t *blockSiteCount;
	uint32_t *blockCompressedLen;
	uint32_t *blockUncompressedLen;

	// for collecting the block weight (number of fluid sites)
	uint32_t *blockWeight;

	// declarations for zlib decompression
	z_stream strm;
	unsigned char *compressedBuffer;
	unsigned char *decompressedBuffer;

	// declarations for site information
	uint32_t siteIsSimulated;
	uint32_t linkType;
	uint32_t siteType;
	float wallDistance;
	uint32_t configEntry;
	uint32_t hasWallNormal;
	float xWallNormal, yWallNormal, zWallNormal;

	// allocations to XDR_DECODE from file stream
	gmyFile = fopen(fname_gmy, "r");
	if (gmyFile == NULL)
	{
		fprintf(stderr, "couldn't open input file %s\n", fname_gmy);
		return 1;
	}
	xdrstdio_create(&xdrs, gmyFile, XDR_DECODE);
	i = 0;

	// allocations for simple zlib decompression
	strm.opaque = Z_NULL;
	strm.zalloc = Z_NULL;
	strm.zfree  = Z_NULL;

	// reading gmy header information
	xdr_u_int(&xdrs, &hlbMagic);
	xdr_u_int(&xdrs, &gmyMagic);
	xdr_u_int(&xdrs, &gmyVersion);
	xdr_u_int(&xdrs, &xBlockCount);
	xdr_u_int(&xdrs, &yBlockCount);
	xdr_u_int(&xdrs, &zBlockCount);
	xdr_u_int(&xdrs, &blockEdgeLen);
	xdr_u_int(&xdrs, &padding);

	DEBUGMSG("    hlbMagic: %u\n", hlbMagic);
	DEBUGMSG("    gmyMagic: %u\n", gmyMagic);
	DEBUGMSG("  gmyVersion: %u\n", gmyVersion);
	DEBUGMSG(" xBlockCount: %u\n", xBlockCount);
	DEBUGMSG(" yBlockCount: %u\n", yBlockCount);
	DEBUGMSG(" zBlockCount: %u\n", zBlockCount);
	DEBUGMSG("blockEdgeLen: %u\n", blockEdgeLen);
	DEBUGMSG("     padding: %u\n", padding);

	// calculating the number of blocks filling the complete volume
	totalBlockCount = xBlockCount * yBlockCount * zBlockCount;
	// calculating the maximum number of sites per block
	maxSitesPerBlock = blockEdgeLen * blockEdgeLen * blockEdgeLen;

	// allocating arrays to read block header information
	blockSiteCount       = new uint32_t[totalBlockCount];
	blockCompressedLen   = new uint32_t[totalBlockCount];
	blockUncompressedLen = new uint32_t[totalBlockCount];
	blockWeight          = new uint32_t[totalBlockCount];

	// reading block header information to arrays
	uint32_t num_entries = 0;
	while (i < totalBlockCount)
	{
		// number of sites actually stored for the block
		xdr_u_int(&xdrs, &blockSiteCount[i]);
		// compressed size of block data in bytes
		xdr_u_int(&xdrs, &blockCompressedLen[i]);
		// uncompressed size of block data in bytes
		xdr_u_int(&xdrs, &blockUncompressedLen[i]);

		// initialise block weights to zero
		blockWeight[i] = 0;

		DEBUGMSG("               block: %u\n", i);
		DEBUGMSG("      blockSiteCount: %u\n", blockSiteCount[i]);
		DEBUGMSG("  blockCompressedLen: %u\n", blockCompressedLen[i]);
		DEBUGMSG("blockUncompressedLen: %u\n", blockUncompressedLen[i]);

		num_entries += blockSiteCount[i];

		// Block data offset can be used to skip directly to blocks / parallel access
		// blockDataOffset[i+1] = blockCompressedLen[i] + blockDataOffset[i];
		i++;
	}

	// reiterate blocks, this time read actual data sets
	for (i = 0; i < totalBlockCount; i++)
	{
		DEBUGMSG("block: %u\n", i);
		if (i % 10000 == 0)
		{
			printf("block %010u of %u\n", i, totalBlockCount);
		}

		// there is only block data available if there is at least one lattice site on the block
		if (blockSiteCount[i] == 0)
		{
			DEBUGMSG("no lattice sites on this block.\n");
		}
		else
		{
			// allocate respective buffers
			compressedBuffer   = new unsigned char[blockCompressedLen[i]];
			decompressedBuffer = new unsigned char[blockUncompressedLen[i]];

			// read compressed block from the file
			ret = fread(compressedBuffer, blockCompressedLen[i], 1, gmyFile);

			// initialise zlib and check if it is ok
			ret = inflateInit(&strm);
			if (ret != Z_OK)
				printf("DECOMPRESSION INIT. ERROR FOR BLOCK\n");

			// specify number of bytes to decompress (inflate in zlib-speech)
			strm.avail_in = blockCompressedLen[i];
			// specify compressed input buffer
			strm.next_in = compressedBuffer;
			// specify number of uncompressed bytes expected after decompression (inflate in zlib-speech)
			// (this is normally used to to read chunks to/from a stream)
			strm.avail_out = blockUncompressedLen[i];
			// specify compressed output buffer
			strm.next_out = decompressedBuffer;

			// run zlib decompression and check if it is ok
			ret = inflate(&strm, Z_FINISH);
			if (ret != Z_STREAM_END)
				printf("DECOMPRESSION ERROR FOR BLOCK\n");

			// finalise/end zlib decompression and check if it is ok
			ret = inflateEnd(&strm);
			if (ret != Z_OK)
				printf("DECOMPRESSION END ERROR FOR BLOCK\n");

			delete(compressedBuffer);

			// create a new XDR_DECODE stream on decompressed buffer
			xdrmem_create(&xdrbs, (char *) decompressedBuffer, blockUncompressedLen[i], XDR_DECODE);

			// loop through all sites on block
			for (uint32_t site = 0; site < maxSitesPerBlock; site++)
			{
				DEBUGMSG("local site index (on block): %u\n", site);

				// check if inside or outside (solid) of the simulation domain
				xdr_u_int(&xdrbs, &siteIsSimulated);
				DEBUGMSG("siteIsSimulated: %u\n", siteIsSimulated);

				if (siteIsSimulated == 1)
				{
					siteType = 0; // default type is FLUID = 0
					uint32_t link = 0;
					uint32_t num_biedges = 0;

					// only consider 26 links (because there is no link to itself)
					for (link = 0; link < 26; link++)
					{
						// check type for every link to our neighbours
						xdr_u_int(&xdrbs, &linkType);
						DEBUGMSG("link %u, linkType: %u\n", link, linkType);
						switch (linkType)
						{
							case 0 :
								// linkType = FLUID
								DEBUGMSG("link %u, FLUID - no further data\n", link);
								break;
							case 1 :
								// linkType = WALL
								// read float of distance from nearest obstacle
								xdr_float(&xdrbs, &wallDistance);
								DEBUGMSG("link %u, WALL - wallDistance: %f\n", link, wallDistance);
								switch (siteType)
								{
									// check for combinations of linkTypes to specify combined siteTypes
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
								// read float of distance from nearest obstacle
								xdr_u_int(&xdrbs, &configEntry);
								// read float of distance from nearest obstacle
								xdr_float(&xdrbs, &wallDistance);
								DEBUGMSG("link %u, INLET - configEntry: %u\n", link, configEntry);
								DEBUGMSG("link %u, INLET - wallDistance: %f\n", link, wallDistance);
								switch (siteType)
								{
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
								// read index of associated configuration in config.xml
								xdr_u_int(&xdrbs, &configEntry);
								// read float of distance from nearest obstacle
								xdr_float(&xdrbs, &wallDistance);
								DEBUGMSG("link %u, OUTLET - configEntry: %u\n", link, configEntry);
								DEBUGMSG("link %u, OUTLET - wallDistance: %f\n", link, wallDistance);
								switch (siteType)
								{
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
								printf("unrecognised linkType %u\n", linkType);
								exit(1);
						}
						if (linkType == 0) // link's a biedge
						{
							num_biedges += 1;
						}
					} // link loop

					// add to block weight according to desired site type weighting
					switch(siteType)
					{
						case 0: // FLUID
							blockWeight[i] += weight_FLUID;
							break;
						case 1: // WALL
							blockWeight[i] += weight_WALL;
							break;
						case 2: // INLET
							blockWeight[i] += weight_INLET;
							break;
						case 3: // OUTLET
							blockWeight[i] += weight_OUTLET;
							break;
						case 4: // INLET+WALL
							blockWeight[i] += weight_INLETWALL;
							break;
						case 5: // OUTLET+WALL
							blockWeight[i] += weight_OUTLETWALL;
							break;
						default:
							printf("bad siteType %u\n", siteType);
							return 1;
					}

					// check if there are wall normal coordinates to be read
					xdr_u_int(&xdrbs, &hasWallNormal);
					DEBUGMSG("hasWallNormal: %u\n", hasWallNormal);
					if (hasWallNormal == 1)
					{
						// read wall normal coordinates as separate floats
						xdr_float(&xdrbs, &xWallNormal);
						xdr_float(&xdrbs, &yWallNormal);
						xdr_float(&xdrbs, &zWallNormal);
						DEBUGMSG("xWallNormal: %f\n", xWallNormal);
						DEBUGMSG("yWallNormal: %f\n", yWallNormal);
						DEBUGMSG("zWallNormal: %f\n", zWallNormal);
					}
					else if (hasWallNormal != 0)
					{
						printf("hasWallNormal missing or wrong... (%u)\n", hasWallNormal);
						exit(1);
					}
				}
				else
				{
					// a solid site does not need a site type
					siteType = -1;
				}
			} // site loop
			// destroy the XDR_DECODE stream on decompressedBuffer
			xdr_destroy(&xdrbs);
			delete(decompressedBuffer);
		} // if block not empty
	}

	// destroy the XDR_DECODE stream on the gmy file
	xdr_destroy(&xdrs);

	// close the gmy file handle
	fclose(gmyFile);

	// write the new gmy+ file
	FILE *gmyplusFile = NULL;
	gmyplusFile = fopen(fname_gmyplus, "w");
	if (gmyplusFile == NULL)
	{
		fprintf(stderr, "couldn't open output file %s\n", fname_gmyplus);
		return 1;
	}
	xdrstdio_create(&xdrsw, gmyplusFile, XDR_ENCODE);

	// write preamble
	xdr_u_int(&xdrsw, &hlbMagic);
	xdr_u_int(&xdrsw, &gmyMagic);
	xdr_u_int(&xdrsw, &gmyVersion);
	xdr_u_int(&xdrsw, &xBlockCount);
	xdr_u_int(&xdrsw, &yBlockCount);
	xdr_u_int(&xdrsw, &zBlockCount);
	xdr_u_int(&xdrsw, &blockEdgeLen);
	xdr_u_int(&xdrsw, &padding);

	// write header with the block weight
	for (i = 0; i < totalBlockCount; i++)
	{
		xdr_u_int(&xdrsw, &blockSiteCount[i]);
		xdr_u_int(&xdrsw, &blockWeight[i]);
		xdr_u_int(&xdrsw, &blockCompressedLen[i]);
		xdr_u_int(&xdrsw, &blockUncompressedLen[i]);
	}

	// reopen input file and skip to the start of the block data section (skip header)
	gmyFile = fopen(fname_gmy, "r");
	if (gmyFile == NULL)
	{
		fprintf(stderr, "couldn't open input file %s\n", fname_gmy);
		return 1;
	}
	xdrstdio_create(&xdrs, gmyFile, XDR_DECODE);
	uint32_t temp;
	uint32_t num_to_skip = 8 + totalBlockCount * 3;
	for (uint32_t skip = 0; skip < num_to_skip; skip++)
	{
		xdr_u_int(&xdrs, &temp);
	}

	// now copy every byte (should make sure we do this more efficiently later)
	printf("copying body of .gmy");
	unsigned char uchar = 0;
	for (;;)
	{
		uchar = fgetc(gmyFile);		// from input
		if (feof(gmyFile))
		{
			break;
		}
		fputc(uchar, gmyplusFile);	// to output
	}

	// close all files
	xdr_destroy(&xdrs);
	xdr_destroy(&xdrsw);
	fclose(gmyFile);
	fclose(gmyplusFile);

	return 0;
}

int main(int argc, char **argv)
{
	if (argc == 9)
	{
		char *fname_gmy = argv[1];
		char *fname_gmyplus = argv[2];
		uint32_t weight_FLUID = atoi(argv[3]);
		uint32_t weight_WALL = atoi(argv[4]);
		uint32_t weight_INLET = atoi(argv[5]);
		uint32_t weight_OUTLET = atoi(argv[6]);
		uint32_t weight_INLETWALL = atoi(argv[7]);
		uint32_t weight_OUTLETWALL = atoi(argv[8]);

		int ret = gmy2gmyplus(fname_gmy, fname_gmyplus, weight_FLUID, weight_WALL, weight_INLET, weight_OUTLET, weight_INLETWALL, weight_OUTLETWALL);
		if (ret != 0)
		{
			printf("error while converting...\n");
			return 1;
		}
		return 0;
	}
	/* otherwise, report usage */
	else
	{
		fputs("Usage: gmy2gmy+ file.gmy file.gmy+ weight_FLUID weight_WALL weight_INLET weight_OUTLET weight_INLETWALL weight_OUTLETWALL\n", stderr);
		return 1;
	}
}
