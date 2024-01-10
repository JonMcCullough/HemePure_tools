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

class VX_Section
{
	public:
		uint64_t num_records;
		uint64_t num_lints_per_record;
		int64_t *data;
		bool chunking_prepared;
		uint64_t chunk, range_min;
		std::map<uint64_t, std::vector<uint64_t>> *block_map;
		std::map<uint64_t, std::vector<uint64_t>>::iterator block_map_it;
		std::map<uint64_t, std::vector<uint64_t>>::iterator block_map_it_test;

		VX_Section(uint64_t num_records, uint64_t num_lints_per_record)
		{
			this->num_records = num_records;
			this->num_lints_per_record = num_lints_per_record;
			data = new int64_t[num_records * num_lints_per_record];
			block_map = new std::map<uint64_t, std::vector<uint64_t>>();
			chunk = 0; range_min = 0;
			this->chunking_prepared = false;
		}

		~VX_Section()
		{
			delete[] data;
		}

		//// read in the pre-specified number of records from the input file starting at the given offset
		//void read_at_offset(MPI_File *in, MPI_Offset offset)
		//{
		//	MPI_Status status;
		//	int count = num_records * num_lints_per_record;
		//	MPI_File_read_at(*in, offset, (void*) data, count, MPI_LONG, &status);
		//}

		void read_at_offset(MPI_File *in, MPI_Offset offset)
		{
			MPI_Status status;
			MPI_Datatype massivetype;
			MPI_Type_contiguous(num_lints_per_record, MPI_LONG, &massivetype);
			MPI_Type_commit(&massivetype);
			if(num_records > INT_MAX) {
				printf("num_records is bigger than signed int (in read_at_offset)\n");
				exit(1);
			}
			int count = num_records;
			MPI_File_read_at(*in, offset, (void*) data, count, massivetype, &status);
		}

		void get_record_xyz(uint64_t record, uint64_t *x, uint64_t *y, uint64_t *z)
		{
			uint64_t x64 = data[record * num_lints_per_record + 0];
			uint64_t y64 = data[record * num_lints_per_record + 1];
			uint64_t z64 = data[record * num_lints_per_record + 2];

			//			assert(x64 <= UINT_MAX);
			//			assert(y64 <= UINT_MAX);
			//			assert(z64 <= UINT_MAX);

			*x = (uint64_t)x64;
			*y = (uint64_t)y64;
			*z = (uint64_t)z64;
		}

		void calc_max_min_blocks(uint64_t blockEdgeLen, uint64_t *xBlockCount, uint64_t *yBlockCount, uint64_t *zBlockCount)
		{
			// calculate the block count
			*xBlockCount = 0;
			*yBlockCount = 0;
			*zBlockCount = 0;
			for(uint64_t i = 0; i < num_records; i++) {

				// get the x, y, z position of this record's site
				uint64_t x, y, z;
				get_record_xyz(i, &x, &y, &z);

				// bin it into its block
				uint64_t blockX = x / blockEdgeLen + 1;
				uint64_t blockY = y / blockEdgeLen + 1;
				uint64_t blockZ = z / blockEdgeLen + 1;

				// get the (local) maximum block grid dimensions
				if(blockX > *xBlockCount) {
					*xBlockCount = blockX;
				}
				if(blockY > *yBlockCount) {
					*yBlockCount = blockY;
				}
				if(blockZ > *zBlockCount) {
					*zBlockCount = blockZ;
				}
				// check for negative x, y or z
				if(x < 0 || y < 0 || z < 0) {
					fprintf(stderr, "ERROR: x, y, z for a site cannot be negative (site %lu has coords %lu,%lu,%lu).\n", i, x, y, z);
					MPI_Finalize();
					exit(1);
				}
			}
		}

		void calc_block_mapping(uint64_t blockEdgeLen, uint64_t xBlockCount, uint64_t yBlockCount, uint64_t zBlockCount)
		{
			for(uint64_t i = 0; i < num_records; i++) {
				// get the x, y, z position of this record's site
				uint64_t x, y, z;
				get_record_xyz(i, &x, &y, &z);

				// bin it into its block
				uint64_t blockX = x / blockEdgeLen;
				uint64_t blockY = y / blockEdgeLen;
				uint64_t blockZ = z / blockEdgeLen;
				uint64_t blockIndex = blockX * yBlockCount * zBlockCount + blockY * zBlockCount + blockZ;

				// check if this block already exists in the map
				auto blockKey = block_map->find(blockIndex);
				if(blockKey == block_map->end()) {
					// if not present, create a vector in which to store record indices
					block_map->emplace(blockIndex, std::vector<uint64_t>());
				}

				// add this record's index to the relevant block
				block_map->at(blockIndex).push_back(i);

				// debug
				if(block_map->at(blockIndex).size() > 512) {
					printf("Block has exceeded 512 sites in calc_block_mapping(), for block %lu [%lu %lu %lu], adding site %lu %lu %lu\n", blockIndex, blockX, blockY, blockZ, x, y, z);
					printf("Current block contents:\n");
					for(uint32_t n = 0; n < block_map->at(blockIndex).size(); n++) {
						uint64_t nx, ny, nz;
						get_record_xyz(block_map->at(blockIndex)[n], &nx, &ny, &nz);
						printf("%lu]] (%lu) %lu %lu %lu\n", blockIndex, block_map->at(blockIndex)[n], nx, ny, nz);
					}
					MPI_Finalize();
					exit(1);
				}
			}
		}

		// print out an ordered list of blocks that (partially or fully) exist on this rank
		void print_blocks_on_rank(int32_t this_rank)
		{
			for(auto it = block_map->begin(); it != block_map->end(); ++it) {
				printf("Rank %d, block %lu\n", this_rank, it->first);
			}
		}

		// this function should be called first, BEFORE any calls to get_next_chunk_of_records()
		void prepare_for_chunking(uint64_t chunk)
		{
			this->block_map_it_test = block_map->begin();
			this->chunking_prepared = true;
			this->chunk = chunk;
			this->range_min = 0;
			this->block_map_it = block_map->begin();
		}

		vector<int64_t> * get_next_chunk_of_records()
		{
			assert(this->chunking_prepared == true);
			std::vector<int64_t> *data_chunk = new std::vector<int64_t>();
			data_chunk->reserve(chunk*159);
			uint64_t range_max = range_min + chunk;
			while(block_map_it != block_map->end()) {
				uint64_t current_block_index = block_map_it->first;
				if(current_block_index >= range_max) {
					break;
				}
				if (current_block_index >= range_min) {
					// add the records in this block to the data_chunk this function returns
					add_records_to_vector(data_chunk, &block_map_it->second);
				} ++block_map_it;
			}
			range_min += chunk;
			return data_chunk;
		}

		uint64_t get_next_chunk_of_records_test()
		{
			assert(this->chunking_prepared == true);
			uint64_t num_records = 0;
			uint64_t range_max = range_min + chunk;
			while(block_map_it_test != block_map->end()) {
				uint64_t current_block_index = block_map_it_test->first;
				if(current_block_index >= range_max) {
					break;
				}
				if (current_block_index >= range_min) {
					num_records += (&block_map_it_test->second)->size();
				} ++block_map_it_test;
			}
			return num_records;
		}

		void add_records_to_vector(std::vector<int64_t> *vec, std::vector<uint64_t> *indices)
		{
			for(auto it = indices->begin(); it != indices->end(); ++it) {
				uint64_t record = *it;
				uint64_t record_start = record * num_lints_per_record;
				uint64_t record_end = record_start + num_lints_per_record;
				// apend the record data to the end of the vector
				vec->insert(vec->end(), &data[record_start], &data[record_end]);
			}
		}

		uint64_t get_num_lints_per_record()
		{
			return num_lints_per_record;
		}
};

typedef struct
{
	uint32_t linkType;
	uint32_t configID;
	float wallDistance;
} VX_LINK;

typedef struct
{
	uint64_t x, y, z;
	bool hasWallNormal;
	float normalX, normalY, normalZ;
	VX_LINK links[26];
} VX_ENTRY;

typedef struct
{
	vector<VX_ENTRY*> data;
} VX;

typedef struct
{
	vector<VX_ENTRY*> *sites;
} BLOCK_ENTRY;

// convert long int to double
double plint_to_double(uint64_t i) {
	//union {double f; uint64_t i;} u;
	//u.i = i;
	//return u.f;
	double result = 0;
	memcpy(&result, &i, sizeof(i));
	return result;
}

// read in the given vx file (binary output from palabos)
VX_Section * read_vx_section(char *fname_vx, MPI_Comm comm)
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

	bool store_records_on_rank0 = false;

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

	// check if the memory needed will be higher than the danger level set in WARN_FOR_MEM_OVER_GB
	double mem_needed_GB = (num_records_on_this_rank * size_of_record)/(1024.0*1024.0*1024.0);
	if(mem_needed_GB > WARN_FOR_MEM_OVER_GB)
		fprintf(stderr, "WARNING: Memory required by rank %05d will be ~%.2fG (WARN LEVEL SET AT %dG).\n", this_rank, mem_needed_GB, WARN_FOR_MEM_OVER_GB);

	// initialise a fresh VX_Section object to contain whatever section of the file we read in
	VX_Section *section = new VX_Section(num_records_on_this_rank, num_lints_per_record);

	// read in the number of records assigned to this rank, starting at the appropriate offset
	DEBUGMSG("Rank %d: reading in %ld/%ld records, starting at offset %lld.\n", this_rank, num_records_on_this_rank, total_num_records, offset);
	section->read_at_offset(&vx_in, offset);
	DEBUGMSG("Rank %d: finished reading. Closing infile.\n", this_rank);
	MPI_File_close(&vx_in);

	return section;
}

VX * convert_records(uint64_t num_lints, int64_t *data)
{
	// create vector to store this rank's entries in
	VX *vx = new VX();

	uint64_t i = 0;
	while(i < num_lints) {

		// initialise a fresh new object to store the next entry
		VX_ENTRY *vx_entry = new VX_ENTRY();
		vx_entry->hasWallNormal = false;
		vx_entry->normalX = 0.0;
		vx_entry->normalY = 0.0;
		vx_entry->normalZ = 0.0;
		uint32_t num_intersections = 0;

		uint64_t x, y, z;
		x = data[i]; i++;
		y = data[i]; i++;
		z = data[i]; i++;

		vx_entry->x = x;
		vx_entry->y = y;
		vx_entry->z = z;
		//		printf("Entry: %u %u %u\n", vx_entry->x, vx_entry->y, vx_entry->z);

		for(uint32_t j = 0; j < 26; j++) {
			uint64_t lu_wallDistance, lu_normalX, lu_normalY, lu_normalZ, linkType, configID;
			float wallDistance=0.0, normalX=0.0, normalY=0.0, normalZ=0.0;

			linkType = data[i]; i++;
			configID = data[i]; i++;
			lu_wallDistance = data[i]; i++;
			lu_normalX = data[i]; i++;
			lu_normalY = data[i]; i++;
			lu_normalZ = data[i]; i++;

			vx_entry->links[j].linkType = (uint32_t) linkType;
			vx_entry->links[j].configID = (uint32_t) configID;

			wallDistance = (float) plint_to_double(lu_wallDistance);
			normalX = (float) plint_to_double(lu_normalX);
			normalY = (float) plint_to_double(lu_normalY);
			normalZ = (float) plint_to_double(lu_normalZ);

			vx_entry->links[j].wallDistance = wallDistance;

			// if link is to a wall...
			if(vx_entry->links[j].linkType == 1) {
				vx_entry->hasWallNormal = true;
				vx_entry->normalX += normalX;
				vx_entry->normalY += normalY;
				vx_entry->normalZ += normalZ;
				num_intersections++;
			}

			//			printf("Link %u = %u %u %f %f %f %f\n", j, vx_entry->links[j].linkType, vx_entry->links[j].configID, vx_entry->links[j].wallDistance, normalX, normalY, normalZ);
		}

		// get average normal
		if(vx_entry->hasWallNormal) {
			vx_entry->normalX /= num_intersections;
			vx_entry->normalY /= num_intersections;
			vx_entry->normalZ /= num_intersections;
		}

		//		printf("hasWallNormal: %d: (%lf, %lf, %lf)\n", vx_entry->hasWallNormal, vx_entry->normalX, vx_entry->normalY, vx_entry->normalZ);

		// store new entry
		vx->data.push_back(vx_entry);
	}
	//	printf("Total number of entries is %ld\n", vx->data.size());

	return vx;
}

void write_chunk_to_file(VX *vx, uint64_t offset, uint64_t num_blocks_in_chunk, uint64_t xBlockCount, uint64_t yBlockCount, uint64_t zBlockCount, uint64_t blockEdgeLen, XDR &xdrs, FILE *gmyFileBlock) {

	XDR xdrbs;
	const uint64_t totalBlockCount = xBlockCount * yBlockCount * zBlockCount;
	const uint64_t maxSitesPerBlock = blockEdgeLen * blockEdgeLen * blockEdgeLen;

	// IMPORTANT: First make sure that num_blocks_in_chunk is reduced if this is the final chunk and it does not add up exactly to totalBlockCount. Otherwise the header will be written with too many blocks in it.
	if(offset + num_blocks_in_chunk >= totalBlockCount) {
		num_blocks_in_chunk = totalBlockCount - offset;
	}

	// initialise the block array for this chunk
	BLOCK_ENTRY *blocks = new BLOCK_ENTRY[num_blocks_in_chunk];
	for(uint64_t i = 0; i < num_blocks_in_chunk; i++) {
		blocks[i].sites = new vector<VX_ENTRY*>();
	}

	// calculate how many sites are in each block
	for(uint64_t i = 0; i < vx->data.size(); i++) {
		uint64_t blockX = vx->data[i]->x / blockEdgeLen;
		uint64_t blockY = vx->data[i]->y / blockEdgeLen;
		uint64_t blockZ = vx->data[i]->z / blockEdgeLen;
		uint64_t blockIndex = blockX * yBlockCount * zBlockCount + blockY * zBlockCount + blockZ;
		if(blockIndex >= totalBlockCount) {
			fprintf(stderr, "ERROR: Block %lu %lu %lu has index %lu, greater than totalBlockCount %lu.\n", blockX, blockY, blockZ, blockIndex, totalBlockCount);
			MPI_Finalize();
			exit(1);
		} else if(blockIndex < offset || blockIndex >= offset + num_blocks_in_chunk) {
			fprintf(stderr, "ERROR: Block %lu %lu %lu has index %lu, outside chunk range [%lu, %lu).\n", blockX, blockY, blockZ, blockIndex, i, i + offset);
			MPI_Finalize();
			exit(1);
		}

		// add site to block (keep its pointer)
		blocks[blockIndex - offset].sites->push_back(vx->data[i]);
	}

	
	for(uint64_t i = 0; i < num_blocks_in_chunk; i++) {
	     // if( blocks[i].sites->size() > 0 ) printf("block %u: blockSiteCount %u\n", i+offset, blocks[i].sites->size());

		// check for problems (e.g. exceeding max number of sites per block)
		if(blocks[i].sites->size() > maxSitesPerBlock) {
			fprintf(stderr, "ERROR: Block %lu exceeds %lu sites (%lu).\n", i + offset, maxSitesPerBlock, blocks[i].sites->size());
			MPI_Finalize();
			exit(1);
		}
	}
	
	// declarations for zlib compression
	z_stream strm;

	uint64_t max_buffer_size = maxSitesPerBlock * (2 + 26 * 4 * sizeof(uint64_t));
	unsigned char *compressedBuffer = new unsigned char [max_buffer_size];
	unsigned char *decompressedBuffer = new unsigned char [max_buffer_size];
	
	// write block data
	VX_ENTRY * blockData[blockEdgeLen][blockEdgeLen][blockEdgeLen];
	for (uint64_t i = 0; i < num_blocks_in_chunk; i++) {
		//		printf("Block %u\n", i);

		// there is only block data available if there is at least one lattice site on the block
		if (blocks[i].sites->size() == 0) {
			// write empty block header entry
			uint32_t blockSiteCount = 0, blockCompressedLen = 0, blockUncompressedLen = 0;
			xdr_u_int(&xdrs, &blockSiteCount);
			xdr_u_int(&xdrs, &blockCompressedLen);
			xdr_u_int(&xdrs, &blockUncompressedLen);

			continue;
		}

		// build block
		for (uint64_t ix = 0; ix < blockEdgeLen; ix++) {
			for (uint64_t iy = 0; iy < blockEdgeLen; iy++) {
				for (uint64_t iz = 0; iz < blockEdgeLen; iz++) {
					blockData[ix][iy][iz] = NULL; // initialise all sites to solid
				}
			}
		}

		// register fluid site pointers within our block
		uint32_t blockUncompressedLen = 0;
		for(vector<VX_ENTRY*>::iterator it = blocks[i].sites->begin(); it != blocks[i].sites->end(); ++it) {
			uint64_t localx = (*it)->x % blockEdgeLen;
			uint64_t localy = (*it)->y % blockEdgeLen;
			uint64_t localz = (*it)->z % blockEdgeLen;
			if(blockData[localx][localy][localz] != NULL) {
				fprintf(stderr, "ERROR: Duplicate site at (local) %lu %lu %lu, for block %lu.\n", localx, localy, localz, i + offset);
				MPI_Finalize();
				exit(1);
			}
			blockData[localx][localy][localz] = *it;
		}

		// work out the uncompressed length of data to be written for this block
		for (uint64_t ix = 0; ix < blockEdgeLen; ix++) {
			for (uint64_t iy = 0; iy < blockEdgeLen; iy++) {
				for (uint64_t iz = 0; iz < blockEdgeLen; iz++) {
					VX_ENTRY *siteptr = blockData[ix][iy][iz];
					if(siteptr == NULL) { // if solid
						blockUncompressedLen += sizeof(uint32_t); // siteIsSimulated
					} else { // if fluid
						blockUncompressedLen += sizeof(uint32_t); // siteIsSimulated
						for(uint32_t m = 0; m < 26; m++) {
							blockUncompressedLen += sizeof(uint32_t); // linkType
							uint32_t linkType = siteptr->links[m].linkType;
							switch (linkType) {
								case 0: // linkType = FLUID (no further data)
									break;
								case 1: // linkType = WALL (write distance to nearest obstacle)
									blockUncompressedLen += sizeof(float); // wallDistance
									break;
								case 2:
								case 3: // linkType = INLET or OUTLET (write config ID and distance to nearest obstacle)
									blockUncompressedLen += sizeof(uint32_t); // configEntry
									blockUncompressedLen += sizeof(float); // wallDistance
									break;
								default:
									fprintf(stderr, "ERROR: Unrecognised linkType %u.\n", linkType);
									MPI_Finalize();
									exit(1);
							}
						}
						uint32_t hasWallNormal = (siteptr->hasWallNormal == true);
						blockUncompressedLen += sizeof(uint32_t); // hasWallNormal
						if (hasWallNormal == 1) {
							blockUncompressedLen += sizeof(float); // normalX
							blockUncompressedLen += sizeof(float); // normalY
							blockUncompressedLen += sizeof(float); // normalZ
						}
					}
				}
			}
		}
		DEBUGMSG("Block: %lu, blockUncompressedLen: %u\n", i + offset, blockUncompressedLen);
		
		if (blockUncompressedLen>65535)
		{
			printf("Block: %lu, blockUncompressedLen: %u\n", i + offset, blockUncompressedLen);
			
			for (uint64_t ix = 0; ix < blockEdgeLen; ix++) {
				for (uint64_t iy = 0; iy < blockEdgeLen; iy++) {
					for (uint64_t iz = 0; iz < blockEdgeLen; iz++) {
						VX_ENTRY *siteptr = blockData[ix][iy][iz];
						if(siteptr == NULL) { // if solid
							printf("Site xyz( %lu, %lu, %lu) is solid", ix, iy, iz); // siteIsSimulated
						} else { // if fluid	
							uint32_t localW = 0;	
							localW += sizeof(uint32_t); // siteIsSimulated
							int w=0;
							int io=0;
							int fl=0;
							for(uint32_t m = 0; m < 26; m++) {
							
								localW += sizeof(uint32_t); // siteIsSimulated
								uint32_t linkType = siteptr->links[m].linkType;
								switch (linkType) {
									case 0: // linkType = FLUID (no further data)
										fl++;
										break;
									case 1: // linkType = WALL (write distance to nearest obstacle)
										localW += sizeof(float); // siteIsSimulated
										w++;	
										break;
									case 2:
									case 3: // linkType = INLET or OUTLET (write config ID and distance to nearest obstacle)
										io++;	
										localW += sizeof(uint32_t); // siteIsSimulated
										localW += sizeof(float); // siteIsSimulated
										break;
									default:
										fprintf(stderr, "ERROR: Unrecognised linkType %u.\n", linkType);
										MPI_Finalize();
										exit(1);
								}
							}
							uint32_t hasWallNormal = (siteptr->hasWallNormal == true);
							localW += sizeof(uint32_t); // siteIsSimulated
							localW += 3*hasWallNormal*sizeof(float); // siteIsSimulated
						printf("Site xyz( %lu, %lu, %lu) is fluid with w, io, fl, wn (%lu, %lu, %lu, %lu) and local block len %lu \n", ix, iy, iz, w, io, fl, hasWallNormal, localW); // siteIsSimulated
						}
					}
				}
			}

		}

		// create a new XDR_ENCODE stream on decompressed buffer (of size blockUncompressedLen)
		xdrmem_create(&xdrbs, (char *) decompressedBuffer, blockUncompressedLen, XDR_ENCODE);

		// loop through all sites on block
		uint64_t num_fluid_sites_written = 0, localSiteID = -1;
		for (uint64_t ix = 0; ix < blockEdgeLen; ix++) {
			for (uint64_t iy = 0; iy < blockEdgeLen; iy++) {
				for (uint64_t iz = 0; iz < blockEdgeLen; iz++) {
					localSiteID++;
					VX_ENTRY *siteptr = blockData[ix][iy][iz];

					// state if inside or outside (aka solid) of the simulation domain (1 or 0)
					uint32_t siteIsSimulated = (siteptr != NULL);
					xdr_u_int(&xdrbs, &siteIsSimulated);

					// write nothing more for solid sites
					if (siteIsSimulated == 0) {
						continue;
					}
					//printf("Block %u, local site id %u (global %u %u %u)\n", i, localSiteID, siteptr->x, siteptr->y, siteptr->z);
					num_fluid_sites_written++;

					for (uint32_t link=0; link < 26; link++) { // only consider 26 links (because there is no link to itself)
						// write type of link
						uint32_t linkType = siteptr->links[link].linkType;
						xdr_u_int(&xdrbs, &linkType);
						switch (linkType) {
							case 0: // linkType = FLUID (no further data)
								break;
							case 1: // linkType = WALL (write distance to nearest obstacle)
								xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
								break;
							case 2: // linkType = INLET (write inletID and distance to nearest obstacle
								xdr_u_int(&xdrbs, &(siteptr->links[link].configID));
								xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
								break;
							case 3: // linkType = OUTLET (write outletID and distance to nearest obstacle
								xdr_u_int(&xdrbs, &(siteptr->links[link].configID));
								xdr_float(&xdrbs, &(siteptr->links[link].wallDistance));
								break;
							default:
								fprintf(stderr, "ERROR: Unrecognised linkType %u.\n", linkType);
								MPI_Finalize();
								exit(1);
						}
					}

					// state if there are wall normal coordinates to be read (1 for yes)
					uint32_t hasWallNormal = (siteptr->hasWallNormal == true);
					xdr_u_int(&xdrbs, &hasWallNormal);
					if (hasWallNormal == 1) {
						// write wall normal coordinates as separate floats
						xdr_float(&xdrbs, &siteptr->normalX);
						xdr_float(&xdrbs, &siteptr->normalY);
						xdr_float(&xdrbs, &siteptr->normalZ);
					}
				}
			}
		} // end of "loop through all sites on block"
		DEBUGMSG("num_fluid_sites_written: %lu\n", num_fluid_sites_written);

		// destroy the XDR_DECODE stream on decompressedBuffer
		xdr_destroy(&xdrbs);

		// setup zlib for compression
		strm.zalloc = Z_NULL;
		strm.zfree = Z_NULL;
		strm.opaque = Z_NULL;

		// input
		strm.avail_in = blockUncompressedLen;
		strm.next_in = decompressedBuffer;

		// output
		strm.avail_out = max_buffer_size;
		strm.next_out = compressedBuffer;

		// perform compression
		uint32_t ret;
		ret = deflateInit(&strm, Z_BEST_COMPRESSION);
		if(ret != Z_OK) {
			fprintf(stderr, "ERROR: zlib deflation init.\n");
			MPI_Finalize();
			exit(1);
		}
		ret = deflate(&strm, Z_FINISH);
		if (ret != Z_STREAM_END) {
			fprintf(stderr, "ERROR: Deflation error for block.\n");
			MPI_Finalize();
			exit(1);
		}
		ret = deflateEnd(&strm);
		if (ret != Z_OK) {
			fprintf(stderr, "ERROR: Deflation end error for block.\n");
			MPI_Finalize();
			exit(1);
		}

		// get new compressed size
		uint32_t blockCompressedLen = (unsigned char*)strm.next_out - (unsigned char*)compressedBuffer;
		DEBUGMSG("Block: %lu, blockUncompressedLen: %u, blockCompressedLen: %u\n", i + offset, blockUncompressedLen, blockCompressedLen);

		// write block header entry
		uint32_t blockSiteCount = (uint32_t)blocks[i].sites->size();
		xdr_u_int(&xdrs, &blockSiteCount);
		xdr_u_int(&xdrs, &blockCompressedLen);
		xdr_u_int(&xdrs, &blockUncompressedLen);

		// write compressed block data to (other) file
		//		uint64_t before = ftell(gmyFileBlock);
		ret = fwrite(compressedBuffer, sizeof(char), blockCompressedLen, gmyFileBlock);
		if(ret != blockCompressedLen) {
			fprintf(stderr, "ERROR: writing compressed block data.\n");
			MPI_Finalize();
			exit(1);
		}
		//		uint64_t after = ftell(gmyFileBlock);
		//		DEBUGMSG("before: %lu after: %lu\n", before, after);
	}

	delete[] decompressedBuffer;
	delete[] compressedBuffer;

	for(uint64_t i = 0; i < num_blocks_in_chunk; i++)
		delete blocks[i].sites;
	delete[] blocks;
}

uint64_t memu()
{
	// get the the current process' status file from the proc filesystem
	FILE* procfile = fopen("/proc/self/status", "r");
	if(procfile == NULL) {
		fprintf(stderr, "Could not open /proc/self/status for reading, in memu().\n");
		MPI_Finalize();
		exit(1);
	}

	long to_read = 8192; char buffer[to_read];
	int read = fread(buffer, sizeof(char), to_read, procfile); read += 1;
	fclose(procfile);

	bool found_vmrss = false;
	char* search_result;
	uint64_t vmrss_kb;

	// look through proc status contents line by line
	char delims[] = "\n";
	char* line = strtok(buffer, delims);
	while (line != NULL && !found_vmrss) {
		search_result = strstr(line, "VmRSS:");
		if (search_result != NULL) {
			sscanf(line, "%*s %ld", &vmrss_kb);
			found_vmrss = true;
		} line = strtok(NULL, delims);
	}
	return vmrss_kb;
}

int vx2gmy(char *fname_vx, char *fname_gmy, uint64_t chunk_of_blocks, MPI_Comm comm)
{
	// GMY header info
	uint32_t hlbMagic=0x686c6221, gmyMagic=0x676d7904, gmyVersion=4, padding=0;
	uint64_t blockEdgeLen = 8;

	// number of blocks (to be determined)
	uint64_t local_xBlockCount=0, local_yBlockCount=0, local_zBlockCount=0, xBlockCount=0, yBlockCount=0, zBlockCount=0;
	uint64_t totalBlockCount=0;

	// read in the voxeliser file data
	VX_Section *section = read_vx_section(fname_vx, comm);

	// find the max and min block counts from the section of the file this rank read in
	section->calc_max_min_blocks(blockEdgeLen, &local_xBlockCount, &local_yBlockCount, &local_zBlockCount);

	int32_t this_rank, num_ranks;
	MPI_Comm_rank(comm, &this_rank);
	MPI_Comm_size(comm, &num_ranks);

	DEBUGMSG("Rank %d: local block counts %lu %lu %lu\n", this_rank, local_xBlockCount, local_yBlockCount, local_zBlockCount);

	// now perform a reduction across all ranks to find the true max block counts in each dimension
	uint64_t local_block_dims[3] = {local_xBlockCount, local_yBlockCount, local_zBlockCount}, global_block_dims[3];
	MPI_Allreduce(local_block_dims, global_block_dims, 3, MPI_UNSIGNED_LONG, MPI_MAX, comm);
	xBlockCount = global_block_dims[0];
	yBlockCount = global_block_dims[1];
	zBlockCount = global_block_dims[2];
	totalBlockCount = xBlockCount * yBlockCount * zBlockCount;

	if(this_rank == 0)
		printf(
				"Rank %d:      global block counts = %lu %lu %lu\n"
				"Rank %d:   total number of blocks = %lu\n", this_rank, xBlockCount, yBlockCount, zBlockCount, this_rank, totalBlockCount);

	// open gmy file for writing with XDR
	FILE *gmyFile = NULL, *gmyFileBlock = NULL;
	XDR xdrs;
	if((gmyFile = fopen(fname_gmy, "w")) == NULL) {
		fprintf(stderr, "Could not open '%s' for writing.\n", fname_gmy);
		exit(1);
	}
	xdrstdio_create(&xdrs, gmyFile, XDR_ENCODE);

	// open temp file for storing compressed block data (since we need to write the header first...)
	string as = string(fname_gmy).append("_blockdata_temp_");
	char *fname_gmy_block = new char[as.length()+1];
	std::strcpy (fname_gmy_block, as.c_str());
	if((gmyFileBlock = fopen(fname_gmy_block, "w")) == NULL) {
		fprintf(stderr, "Could not open '%s' for writing.\n", fname_gmy_block);
		exit(1);
	}

	// Convert 64->32 for gmy format
	uint32_t xBlockCount32 = (uint32_t)xBlockCount;
	uint32_t yBlockCount32 = (uint32_t)yBlockCount;
	uint32_t zBlockCount32 = (uint32_t)zBlockCount;
	uint32_t blockEdgeLen32 = (uint32_t)blockEdgeLen;

	// write gmy header information
	xdr_u_int(&xdrs, &hlbMagic);
	xdr_u_int(&xdrs, &gmyMagic);
	xdr_u_int(&xdrs, &gmyVersion);
	xdr_u_int(&xdrs, &xBlockCount32);
	xdr_u_int(&xdrs, &yBlockCount32);
	xdr_u_int(&xdrs, &zBlockCount32);
	xdr_u_int(&xdrs, &blockEdgeLen32);
	xdr_u_int(&xdrs, &padding);

	// bin all the sites this rank knows about into blocks
	section->calc_block_mapping(blockEdgeLen, xBlockCount, yBlockCount, zBlockCount);

	//section->print_blocks_on_rank(this_rank);

	// To keep memory usage low, rank 0 now has to fetch blocks in "chunks" i.e. by collecting records in the first 100 blocks from all ranks, then the next 100 blocks, and so on, writing the gmy to disk as it goes.
	section->prepare_for_chunking(chunk_of_blocks);
	for(uint64_t i = 0; i < totalBlockCount; i += chunk_of_blocks) {
		if(this_rank == 0)
			printf("Rank 0: want blocks with indices in the range %09lu to %09lu.\n", i, i + chunk_of_blocks);

		uint64_t mmryusge = memu();
		uint64_t num_recd = section->get_next_chunk_of_records_test();

//		mmryusge += num_recd*section->get_num_lints_per_record()*8/1000;
//		fprintf(stderr, "WARNING: Memory required by rank %05d is %lu kB, with num_recd = %lu.\n", this_rank, mmryusge, num_recd);

		// prepare vector of relevant records to send
		std::vector<int64_t> *records = section->get_next_chunk_of_records();

		//printf("Rank %d: has %lu sites.\n", this_rank, records->size());
		//for(uint64_t i=0; i < records->size(); i+=159) {
		//	printf("[%ld %ld %ld] ", (*records)[i], (*records)[i+1], (*records)[i+2]);
		//}
		//printf("\n");

		// get receive counts from each rank
		int *recvcounts = NULL;
		if(this_rank == 0) {
			recvcounts = new int[num_ranks];
		}
		int num_lints_to_send = records->size();
		MPI_Gather(&num_lints_to_send, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, comm);

		// prepare the arrays rank 0 will need to receive all other ranks' records (from current chunk)
		int *displs = NULL;
		int64_t *all_records = NULL;
		uint64_t total_num_lints = 0;
		if(this_rank == 0) {
			displs = new int[num_ranks];
			for(int j = 0; j < num_ranks; j++) {
				displs[j] = total_num_lints;
				total_num_lints += recvcounts[j];
			}
			all_records = new int64_t[total_num_lints];
		}

//		if(this_rank == 0) {
//			mmryusge = memu();
//			fprintf(stderr, "WARNING: Memory required by rank %05d is %lu kB, with total_num_lints = %lu.\n", this_rank, mmryusge, total_num_lints);
//		}

		int64_t *to_send = NULL;
		if(num_lints_to_send > 0) {
			to_send = &(*records)[0]; // this looks so awful...
		}

		//if(this_rank == 0) {
		//	DEBUGMSG("%u] Require %lu records. Break-down:\n", i, total_num_lints/section->get_num_lints_per_record());
		//	for(int j = 0; j < num_ranks; j++) {
		//		DEBUGMSG("%u] Rank %d will send %ld\n", i, j, recvcounts[j]/section->get_num_lints_per_record());
		//	}
		//}

		// collect all relevant records (in this chunk) from all ranks
		MPI_Gatherv(to_send, num_lints_to_send, MPI_LONG, all_records, recvcounts, displs, MPI_LONG, 0, comm);

		// convert, compress and write to file this chunk of blocks
		VX *vx = NULL;
		if(this_rank == 0) {
			vx = convert_records(total_num_lints, all_records);
			write_chunk_to_file(vx, i, chunk_of_blocks, xBlockCount, yBlockCount, zBlockCount, blockEdgeLen, xdrs, gmyFileBlock);

			for(uint64_t j = 0; j < vx->data.size(); j++)
				delete vx->data[j];
		}

		// clean up after this chunk
		delete vx;
		delete[] recvcounts;
		delete[] all_records;
		delete records;

//		if(this_rank == 0) {
//			mmryusge = memu();
//			fprintf(stderr, "WARNING: Memory required by rank %05d is %lu kB.\n", this_rank, mmryusge);
//		}
	}
	xdr_destroy(&xdrs);
	fclose(gmyFile);
	fclose(gmyFileBlock);

	// all done, no problems
	return 0;
}

int main(int argc, char **argv)
{
	assert (sizeof(uint64_t) == sizeof(double));

	if (argc != 3 && argc != 4) {
		fprintf(stderr, "Usage: vx2gmy file.vx file.gmy [BLOCKCHUNKSIZE].\n");
		return 1;
	}

	char *fname_vx = argv[1];
	char *fname_gmy = argv[2];
	uint64_t chunk_of_blocks = 10;
	if(argc == 4) {
		chunk_of_blocks = atoi(argv[3]);
	}

	// set up MPI
	MPI_Init(NULL, NULL);

	uint32_t ret = vx2gmy(fname_vx, fname_gmy, chunk_of_blocks, MPI_COMM_WORLD);
	if (ret != 0) {
		fprintf(stderr, "Error ocurred.\n");
		return 1;
	}

	MPI_Finalize();

	return 0;
}
