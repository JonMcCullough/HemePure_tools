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
#include <cstdio> 

#define DEBUG 0

#if DEBUG==1
#define DEBUGMSG(...) printf(__VA_ARGS__)
#else
#define DEBUGMSG(...)
#endif

// convert long int to double
double plint_to_double(uint64_t i) {
    double result = 0;
    memcpy(&result, &i, sizeof(i));
    return result;
}

void record_to_ascii(int64_t *data) {
    uint32_t i = 0;
    
    uint64_t x, y, z;

    x = data[i]; i++;
    y = data[i]; i++;
    z = data[i]; i++;

    printf("x y z = %lld %lld %lld\n", x, y, z);

    for(uint32_t j = 0; j < 26; j++) {
        uint64_t lu_wallDistance, lu_normalX, lu_normalY, lu_normalZ;
        uint32_t linkType, configID;
        float wallDistance=0.0, normalX=0.0, normalY=0.0, normalZ=0.0;

        linkType = (uint32_t) data[i]; i++;
        configID = (uint32_t) data[i]; i++;
        lu_wallDistance = data[i]; i++;
        lu_normalX = data[i]; i++;
        lu_normalY = data[i]; i++;
        lu_normalZ = data[i]; i++;

        wallDistance = (float) plint_to_double(lu_wallDistance);
        normalX = (float) plint_to_double(lu_normalX);
        normalY = (float) plint_to_double(lu_normalY);
        normalZ = (float) plint_to_double(lu_normalZ);

        printf("Link %u: linkType=%u configID=%u wallDist=%f normal=(%f %f %f)\n", j, linkType, configID, wallDistance, normalX, normalY, normalZ);
    }
}




int main(int argc, char **argv)
{
	assert (sizeof(uint64_t) == sizeof(double));

	if (argc != 2) {
		fprintf(stderr, "Usage: vx2ascii file.vx\n");
		return 1;
	}

	char *fname_vx = argv[1];

    const int64_t num_lints_per_record = 3 + 26 * 6;
    const int64_t size_of_record = num_lints_per_record * sizeof(uint64_t);
    int64_t record[size_of_record];

    FILE *vxfile = fopen(fname_vx, "r");
    if(vxfile == NULL) {
        fprintf(stderr, "Could not open file for reading.\n");
        return 1;
    }

    uint64_t i = 0;
    while( fread(record, 1, size_of_record, vxfile) == size_of_record ) {
        printf("Site %lld :\n", i);
        record_to_ascii(record);
        printf("\n");
        i++;
    }

    fclose(vxfile);

	return 0;
}
