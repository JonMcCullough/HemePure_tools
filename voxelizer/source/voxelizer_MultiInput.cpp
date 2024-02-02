#include "palabos3D.h"
#include "palabos3D.hh"

#include <cmath>
#include <cstdio>
#include <limits>
#include <unordered_map>

#include <iostream>
#include <fstream>

#include "utils.h"

// This is for timing functions -- will check to see if Palabos has them
#include <mpi.h>

// ========== DO NOT TOUCH (START) ========== //
//  BJ: One can actually touch it but not recommended -- this code is needed for the parallel IO 
//  the key functionality here is the send, which trawls all the internal blocks of the 
// multi-block structure and collects the data fields for transfer
namespace plb {
	template<typename T, int nDim>
		class SparseTensorFieldDataTransfer3D : public BlockDataTransfer3D {
			public:
				SparseTensorFieldDataTransfer3D();
				virtual void setBlock(AtomicBlock3D& block);
				virtual void setConstBlock(AtomicBlock3D const& block);
				virtual SparseTensorFieldDataTransfer3D<T,nDim>* clone() const;
				virtual plint staticCellSize() const;
				// Send data from the block into a byte-stream.
				virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
				// Receive data from a byte-stream into the block.
				virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind);
				virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset) {
					receive(domain, buffer, kind);
				}
				virtual void receive(Box3D domain, std::vector<char> const& buffer,
						modif::ModifT kind, std::map<int,std::string> const& foreignIds)
				{
					receive(domain, buffer, kind);
				}
				// Attribute data between two blocks.
				virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
						AtomicBlock3D const& from, modif::ModifT kind);
				virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
						AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset)
				{
					attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
				}
			private:
				TensorField3D<T,nDim>* field;
				TensorField3D<T,nDim> const* constField;
		};

	// Class SparseTensorFieldDataTransfer3D

	template<typename T, int nDim>
		SparseTensorFieldDataTransfer3D<T,nDim>::SparseTensorFieldDataTransfer3D()
		: field(0),
		constField(0)
	{ }

	template<typename T, int nDim>
		void SparseTensorFieldDataTransfer3D<T,nDim>::setBlock(AtomicBlock3D& block) {
			field = dynamic_cast<TensorField3D<T,nDim>*>(&block);
			PLB_ASSERT(field);
			constField = field;
		}

	template<typename T, int nDim>
		void SparseTensorFieldDataTransfer3D<T,nDim>::setConstBlock(AtomicBlock3D const& block) {
			constField = dynamic_cast<TensorField3D<T,nDim> const*>(&block);
			PLB_ASSERT(constField);
		}

	template<typename T, int nDim>
		SparseTensorFieldDataTransfer3D<T,nDim>* SparseTensorFieldDataTransfer3D<T,nDim>::clone() const {
			return new SparseTensorFieldDataTransfer3D<T,nDim>(*this);
		}

	template<typename T, int nDim>
		plint SparseTensorFieldDataTransfer3D<T,nDim>::staticCellSize() const {
			return nDim*sizeof(T);
		}

	template<typename T, int nDim>
		void SparseTensorFieldDataTransfer3D<T,nDim>::send(
				Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const
		{
			PLB_PRECONDITION(constField);
			PLB_PRECONDITION(contained(domain, constField->getBoundingBox()));
			PLB_ASSERT(nDim>=1);

			plint cellSize = staticCellSize();
			buffer.clear();
			std::vector<char> tmpData(cellSize);

			for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
				for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
					for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
						if (constField->get(iX,iY,iZ)[0] != -1) {
							memcpy((void*)(&tmpData[0]), (const void*)(&constField->get(iX,iY,iZ)[0]), nDim*sizeof(T));
							buffer.insert(buffer.end(), tmpData.begin(), tmpData.end());
						}
					}
				}
			}
		}

	template<typename T, int nDim>
		void SparseTensorFieldDataTransfer3D<T,nDim>::receive(
				Box3D domain, std::vector<char> const& buffer, modif::ModifT kind)
		{
			PLB_PRECONDITION(field);
			PLB_PRECONDITION(contained(domain, field->getBoundingBox()));
			PLB_PRECONDITION(domain.nCells()*staticCellSize() == (plint)buffer.size());

			// Avoid dereferencing uninitialized pointer.
			if (buffer.empty()) return;
			plint iData=0;
			for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
				for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
					for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
						memcpy((void*)(&field->get(iX,iY,iZ)[0]), (const void*)(&buffer[iData]), nDim*sizeof(T));
						iData += nDim*sizeof(T);
					}
				}
			}
		}

	template<typename T, int nDim>
		void SparseTensorFieldDataTransfer3D<T,nDim>::attribute(
				Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
				AtomicBlock3D const& from, modif::ModifT kind)
		{
			PLB_PRECONDITION(typeid(from) == typeid(TensorField3D<T,nDim> const&));
			PLB_PRECONDITION(contained(toDomain, field->getBoundingBox()));

			TensorField3D<T,nDim> const& fromField = (TensorField3D<T,nDim> const&) from;
			for (plint iX=toDomain.x0; iX<=toDomain.x1; ++iX) {
				for (plint iY=toDomain.y0; iY<=toDomain.y1; ++iY) {
					for (plint iZ=toDomain.z0; iZ<=toDomain.z1; ++iZ) {
						for (int iDim=0; iDim<nDim; ++iDim) {
							field->get(iX,iY,iZ)[iDim] = fromField.get(iX+deltaX,iY+deltaY,iZ+deltaZ)[iDim];
						}
					}
				}
			}
		}
}

// ========== DO NOT TOUCH (END) ========== //

using namespace plb;

#define DESCRIPTOR descriptors::D3Q19Descriptor


const plint blockSize = 32; 			// Zero means: no sparse representation.
const plint extendedEnvelopeWidth = 2;	// Because the Guo off lattice boundary condition needs 2-cell neighbor access.
const plint borderWidth = 1;			// Because the Guo boundary condition acts in a one-cell layer.
										//		Requirement: margin >= borderWidth.
plint extraLayer = 0;					// Make the bounding box larger; for visualization purposes only.
										//		For the simulation, it is OK to have extraLayer = 0.
const plint margin = 5;					// Extra margin of allocated cells around the obstacle.

bool linkFileUseInnerBorder = true;
bool linkFileUseOuterBorder = false;

// typedefs
typedef double T;
typedef TriangleSet<T>::Triangle Triangle;

// structure which defines an opening
// ----------------------------------
// the surface geometry (*.stl) contains openings, which represent inlets and outlets
template<typename T>
struct Opening {
	Array<T,3> center;
	Array<T,3> normal;
	T innerRadius;
	T outerRadius;
	T area; //JM
	int type; // type = 2 is inlet; type = 3 is outlet
	int index;
}; std::vector<Opening<T> > openings;
// ---

// input variables from .xml
// -------------------------
pluint numIlets, numOlets;
plint referenceDirection = 0;
T dxREL = 0;
T dxABS = 0;

std::vector<Array<T,3> > analysisPoints;
// ---

// variables relating to surface
// -----------------------------
TriangleSet<T>* arteryTriangleSet = 0;
TriangularSurfaceMesh<T>* arterySurfaceMesh = 0;

std::string arteryMeshFileName;
std::unordered_map<plint,int> lidTriangles;
// ---

// ensure we can safely convert double to plint
static_assert(sizeof(T) == sizeof(plint), "ERROR: sizeof(double) != sizeof(plint)!");

// read the user input XML file provided at the command-line
void readGeneralParameters(XMLreader const& generalXML)
{
	generalXML["referenceDirection"].read(referenceDirection);
	generalXML["resolution"].read(dxREL);
	generalXML["DX"].read(dxABS);
	generalXML["stl"].read(arteryMeshFileName);
	generalXML["analysisPoints"]["numIlets"].read(numIlets);
	generalXML["analysisPoints"]["numOlets"].read(numOlets);
	for (pluint i = 1; i <= numIlets+numOlets; ++i) {
		std::vector<T> nextPoint;
		generalXML["analysisPoints"]["point"][i].read(nextPoint);

		if (nextPoint.size() != 3) {
			throw PlbIOException("-> analysis points must be 3-component coordinates");
		}
		analysisPoints.push_back(Array<T,3>(nextPoint[0],nextPoint[1],nextPoint[2]));
	}

	// the surface geometry (stl file) is read into data structures comprised of a set of triangles
	// the DBL constant means that double precision accuracy will be used
	arteryTriangleSet = new TriangleSet<T>(arteryMeshFileName, DBL);
}

/* BJ:  I refactored this to be a general BoxProcessingFunctional 
  so that I can pass it an extra argument besides the links and the flags.
  I also want to pass in a pre-created triangle hash. The only complication is
  that instead of a process() function we now neeed to use processGenericBlocks()
  instead. This takes a vector of multi-block structures as its argument. In our 
  case all 3 multi-block arguments should have an identical multi-block structure
  -- that of the voxel matrix 
  */


template<typename T1, typename T2>
class VoxelToLink3D : public BoxProcessingFunctional3D {
public:
		VoxelToLink3D(Box3D bbox);

		virtual void processGenericBlocks( Box3D domain, std::vector<AtomicBlock3D*> fields);

		virtual VoxelToLink3D<T1,T2>* clone() const {
			return new VoxelToLink3D(*this);
		}
		virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
			modified[0] = modif::nothing;
			modified[1] = modif::staticVariables;
		}
	private:
		plint nX, nY, nZ;
};

// class VoxelToLink3D
// -------------------

template<typename T1, typename T2>
VoxelToLink3D<T1,T2>::VoxelToLink3D(Box3D bbox):
	nX(bbox.x1+1),
	nY(bbox.y1+1),
	nZ(bbox.z1+1)
{ }

// This is the main processor
template<typename T1, typename T2>
void VoxelToLink3D<T1,T2>::processGenericBlocks (Box3D domain, std::vector<AtomicBlock3D*> fields)
//		Box3D domain, ScalarField3D<T1>& flags, TensorField3D<T2,159>& links)
{
	// Unwrap the fields, dereference the pointers and dynamic cast. 
	ScalarField3D<T1>& flags = dynamic_cast<ScalarField3D<T1>&>(*fields[0]);
	TensorField3D<T2,159>& links = dynamic_cast<TensorField3D<T2,159>&>(*fields[1]);

	// This is a hash container of nearby triangles computed the same way as for voxelization
 	AtomicContainerBlock3D& triHashContainer = dynamic_cast<AtomicContainerBlock3D&>(*fields[2]);

	Dot3D absPos = flags.getLocation();
	Dot3D ofs = computeRelativeDisplacement(flags, links);
	
	Array<T ,3> displacement;
	Array<T2,2> xRange, yRange, zRange;
	Array<T2,3> point1, point2;

	T          distce, tmpDistce;
	Array<T,3> intrsc, tmpIntrsc;
	Array<T,3> normal, tmpNormal;

	normal[0] = -2.0;
	normal[1] = -2.0;
	normal[2] = -2.0;

	pluint iWrite, opening, result, whichT;

	//// check conversion of double to plint...
	//distce = -atan(1.0)*4.0;
	//memcpy(&result, &distce, sizeof(distce));
	//distce = -1.0;
	////... and back
	//memcpy(&distce, &result, sizeof(result));
	//if (distce != -atan(1.0)*4.0)
	//	std::runtime_error("--> did not successfully convert to double to plint!");

	distce = std::numeric_limits<T>::max();

	// loop through all sites (inside and outside)
	for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
		for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
			for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
				Array<T2,159>& linkdata = links.get(iX+ofs.x, iY+ofs.y, iZ+ofs.z);

				// site has been marked (inside)
				if (flags.get(iX,iY,iZ)) {
					linkdata[0] = iX+absPos.x;
					linkdata[1] = iY+absPos.y;
					linkdata[2] = iZ+absPos.z;
					iWrite = 2;

					// loop through neighbours of site
					for (plint dx = -1; dx <= +1; ++dx)
						for (plint dy = -1; dy <= +1; ++dy)
							for (plint dz = -1; dz <= +1; ++dz) {
								if (!(dx == 0 && dy == 0 && dz == 0)) {

									if (flags.get(iX+dx,iY+dy,iZ+dz)) iWrite += 6;	// neighbour is a fluid site
									else {											// neighbour is a solid site

										// construct segment (line between p1 and p2)
										point1[0] = linkdata[0]; point2[0] = linkdata[0]+dx;
										point1[1] = linkdata[1]; point2[1] = linkdata[1]+dy;
										point1[2] = linkdata[2]; point2[2] = linkdata[2]+dz;

										displacement[0] = dx;
										displacement[1] = dy;
										displacement[2] = dz;

										// Now we will only get the triangles 
										// that are near us. We use the triangleHasn here
										// 
										Array<T2,2> xRange(
												std::min(point1[0], point2[0]),
												std::max(point1[0], point2[0]));
										Array<T2,2> yRange(
												std::min(point1[1], point2[1]),
												std::max(point1[1], point2[1]));
										Array<T2,2> zRange(
												std::min(point1[2], point2[2]),
												std::max(point1[2], point2[2]));

										TriangleHash<T> triangleHash(triHashContainer);
    									std::vector<plint> possibleTriangles;
    									triangleHash.getTriangles(xRange, yRange, zRange, possibleTriangles);

										std::vector<T> crossings;
										for (pluint iPossible = 0; iPossible < possibleTriangles.size(); ++iPossible) {
											plint iTriangle = possibleTriangles[iPossible];
											if (arterySurfaceMesh->pointOnTriangle(
														point1, point2, 0, iTriangle,
														tmpIntrsc, tmpNormal, tmpDistce) == 1) {
												if (!util::fpequal_abs(tmpDistce, T(), arterySurfaceMesh->eps1)) {
													crossings.push_back(tmpDistce);
												}
												if (crossings.size() == 1 || tmpDistce < distce) {
													distce = tmpDistce; intrsc = tmpIntrsc; normal = tmpNormal;
													whichT = iTriangle;
												}
											}
										}

										bool hasCrossed = true;
										if (crossings.size() == 0)
											std::runtime_error("--> could not find intersection with any triangle (1)!");
										else {
											for (pluint iCrossing = 1; iCrossing < crossings.size(); ++iCrossing) {
												const T eps1 = std::numeric_limits<T>::epsilon()*1.e4;
												if (std::fabs(crossings[iCrossing]-crossings[iCrossing-1]) > eps1)
													hasCrossed = !hasCrossed;
											}
										}
										if (!hasCrossed)
											std::runtime_error("--> could not find intersection with any triangle (2)!");

										// distance should be a fraction of lattice vector
										distce /= norm(displacement);

										T normD = norm(displacement); displacement /= normD;
										T denom = dot(normal, displacement);
										if (denom >= 0.0)
											std::runtime_error("--> outward-pointing normal!");

										if (lidTriangles.find(whichT) != lidTriangles.end()) {
											//std::cout << "---> belonging to lid: " <<
											//lidTriangles.at(whichTriangle) << std::endl;
											linkdata[++iWrite] = openings[
												lidTriangles.at(whichT)].type;
								
											opening = openings[
												lidTriangles.at(whichT)].index;

											
											if (opening < numIlets) {
										
												linkdata[++iWrite] = opening;
												if (opening < 0) throw std::runtime_error("--> invalid opening index (1)!");
											} else {
												opening -= numIlets;
												
												linkdata[++iWrite] = opening;
												if (opening < 0) throw std::runtime_error("--> invalid opening index (2)!");
											}
										} else { // intersection with wall
											linkdata[++iWrite] = 1; ++iWrite;
										}

										// distance to intersection...
										// ... and three components of normal
										memcpy(&result, &distce,    sizeof(distce));
										linkdata[++iWrite] = result;
										memcpy(&result, &normal[0], sizeof(normal[0]));
										linkdata[++iWrite] = result;
										memcpy(&result, &normal[1], sizeof(normal[1]));
										linkdata[++iWrite] = result;
										memcpy(&result, &normal[2], sizeof(normal[2]));
										linkdata[++iWrite] = result;
									}
								}
							}
				} else linkdata[0] = -1;
			}
		}
	}
}

void memu()
{
	// get the the current process' status file from the proc filesystem
	FILE* procfile = fopen("/proc/self/status", "r");

	long to_read = 8192; char buffer[to_read];
	int read = fread(buffer, sizeof(char), to_read, procfile);
	fclose(procfile);

	bool found_vmrss = false;
	char* search_result;
	long vmrss_kb;

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

	std::stringstream output;
	long vmrss_kb_global;
	MPI_Allreduce(&vmrss_kb, &vmrss_kb_global,1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD );

	pcout << "-| memory use: this node: " << vmrss_kb << " kB.  Global max: " << vmrss_kb_global << " kB\n";
}

void writeXML(const pluint num_openings, double dx, double shift_x, double shift_y, double shift_z)
{
	std::ofstream myfile;
	myfile.open ("iolets_block_inputxml_PP.txt");

	myfile << "  <inlets>\n";

	int inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 2) { // inlets
			inletnumber++;
			myfile <<
				"    <inlet>\n"
				"      <!-- index value=" << inletnumber-1 << "-->\n"
				"      <condition subtype=\"cosine\" type=\"pressure\">\n"
				"        <amplitude units=\"mmHg\" value=\"0.0\"/>\n"
				"        <mean units=\"mmHg\" value=\"0.0\"/>\n"
				"        <phase units=\"rad\" value=\"0.0\"/>\n"
				"        <period units=\"s\" value=\"1\"/>\n"
				"        <radius value=\"" << 0.5*(openings[i].innerRadius + openings[i].outerRadius)*dx << "\" units=\"m\"/>\n"
				"        <area value=\"" << openings[i].area*dx*dx << "\" units=\"m^2\"/>\n"
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </inlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </inlet>\n";
		}
	}
	myfile << "  </inlets>\n";
	myfile << "  <outlets>\n";
	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 3) { // outlets
			inletnumber++;
			myfile <<
				"    <outlet>\n"
				"      <!-- index value=" << inletnumber-1 << " -->\n"
				"      <condition subtype=\"cosine\" type=\"pressure\">\n"
				"        <amplitude units=\"mmHg\" value=\"0.0\"/>\n"
				"        <mean units=\"mmHg\" value=\"0.0\"/>\n"
				"        <phase units=\"rad\" value=\"0.0\"/>\n"
				"        <period units=\"s\" value=\"1\"/>\n"
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </outlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </outlet>\n";
		}
	}
	myfile << "  </outlets>\n";

	myfile.close();

	
	myfile.open ("iolets_block_inputxml_VP.txt");

	myfile << "  <inlets>\n";

	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 2) { // inlets
			inletnumber++;
			myfile <<
				"    <inlet>\n"
				"      <!-- index value=" << inletnumber-1 << "-->\n"
				"      <condition type=\"velocity\" subtype=\"parabolic\">\n"
				"        <radius value=\"" << 0.5*(openings[i].innerRadius + openings[i].outerRadius)*dx << "\" units=\"m\"/>\n"
				"        <area value=\"" << openings[i].area*dx*dx << "\" units=\"m^2\"/>\n"
				"        <maximum value=\"CHANGE\" units=\"m/s\"/>\n"
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </inlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </inlet>\n";
		}
	}
	myfile << "  </inlets>\n";
	myfile << "  <outlets>\n";
	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 3) { // outlets
			inletnumber++;
			myfile <<
				"    <outlet>\n"
				"      <!-- index value=" << inletnumber-1 << " -->\n"
				"      <condition subtype=\"cosine\" type=\"pressure\">\n"
				"        <amplitude units=\"mmHg\" value=\"0.0\"/>\n"
				"        <mean units=\"mmHg\" value=\"0.0\"/>\n"
				"        <phase units=\"rad\" value=\"0.0\"/>\n"
				"        <period units=\"s\" value=\"1\"/>\n"
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </outlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </outlet>\n";
		}
	}
	myfile << "  </outlets>\n";

	myfile.close();
	

	myfile.open ("iolets_block_inputxml_VfP.txt");

	myfile << "  <inlets>\n";

	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 2) { // inlets
			inletnumber++;
			myfile <<
				"    <inlet>\n"
				"      <!-- index value=" << inletnumber-1 << "-->\n"
				"      <condition type=\"velocity\" subtype=\"file\">\n"
				"        <path value=\"MESHX_INLET" << inletnumber-1 << "_VELOCITY.txt\"/>\n"
				"        <radius value=\"" << 0.5*(openings[i].innerRadius + openings[i].outerRadius)*dx << "\" units=\"m\"/>\n"
				"        <area value=\"" << openings[i].area*dx*dx << "\" units=\"m^2\"/>\n"
				"       <!-- <maximum value=\"CHANGE (Coupled only)\" units=\"m/s\"/> -->\n"	
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </inlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </inlet>\n";
		}
	}
	myfile << "  </inlets>\n";
	myfile << "  <outlets>\n";
	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 3) { // outlets
			inletnumber++;
			myfile <<
				"    <outlet>\n"
				"      <!-- index value=" << inletnumber-1 << " -->\n"
				"      <condition subtype=\"cosine\" type=\"pressure\">\n"
				"        <amplitude units=\"mmHg\" value=\"0.0\"/>\n"
				"        <mean units=\"mmHg\" value=\"0.0\"/>\n"
				"        <phase units=\"rad\" value=\"0.0\"/>\n"
				"        <period units=\"s\" value=\"1\"/>\n"
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </outlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </outlet>\n";
		}
	}
	myfile << "  </outlets>\n";

	myfile.close();

	
	myfile.open ("iolets_block_inputxml_VfWKf.txt");

	myfile << "  <inlets>\n";

	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 2) { // inlets
			inletnumber++;
			myfile <<
				"    <inlet>\n"
				"      <!-- index value=" << inletnumber-1 << "-->\n"
				"      <condition type=\"velocity\" subtype=\"file\">\n"
				"        <path value=\"MESHX_INLET" << inletnumber-1 << "_VELOCITY.txt\"/>\n"	
				"        <radius value=\"" << 0.5*(openings[i].innerRadius + openings[i].outerRadius)*dx << "\" units=\"m\"/>\n"
				"        <area value=\"" << openings[i].area*dx*dx << "\" units=\"m^2\"/>\n"
				"       <!-- <maximum value=\"CHANGE (Coupled only)\" units=\"m/s\"/> -->\n"
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </inlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </inlet>\n";
		}
	}
	myfile << "  </inlets>\n";
	myfile << "  <outlets>\n";
	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 3) { // outlets
			inletnumber++;
			myfile <<
				"    <outlet>\n"
				"      <!-- index value=" << inletnumber-1 << " -->\n"
				"      <condition subtype=\"fileWK\" type=\"pressure\">\n"
				"        <path value=\"OUTLET" << inletnumber-1 << "_WK.txt.weights.txt\"/>\n"
				"        <R value=\"CHANGE\" units=\"kg/m^4*s\"/>\n"
				"        <C value=\"CHANGE\" units=\"m^4*s^2/kg\"/>\n"
				"        <radius value=\"" << 0.5*(openings[i].innerRadius + openings[i].outerRadius)*dx << "\" units=\"m\"/>\n"
				"        <area value=\"" << openings[i].area*dx*dx << "\" units=\"m^2\"/>\n"
				"      </condition>\n"
				"      <normal units=\"dimensionless\" value=\"(" <<
				openings[i].normal[0] << "," << openings[i].normal[1] << "," << openings[i].normal[2] << ")\"/>\n"
				// Lattice Units
				"      <position units=\"lattice\" value=\"(" <<
				openings[i].center[0] << "," << openings[i].center[1] << "," << openings[i].center[2] << ")\"/>\n"
				"    </outlet>\n";
				//// Physical Units
				// "      <position units=\"m\" value=\"(" <<
				// (openings[i].center[0] - shift_x)*dx << "," << (openings[i].center[1] - shift_y)*dx << "," << (openings[i].center[2] - shift_z)*dx << ")\"/>\n"
				// "    </outlet>\n";
		}
	}
	myfile << "  </outlets>\n";

	myfile.close();
	pcout << "-> writing XML\n";

	// Radius file (openings has this info)
	myfile.open("inlets_radius.txt");
	//myfile << "inlets\n";

	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 2) { // inlets
			inletnumber++;
			myfile << inletnumber << "," << 0.5*(openings[i].innerRadius + openings[i].outerRadius)*dx << "," << openings[i].area*dx*dx << "\n";
		}
	}

	myfile.close();
	myfile.open("outlets_radius.txt");

	inletnumber = 0;
	for (pluint i = 0; i < num_openings; ++i) {
		if (openings[i].type == 3) { // outlets
			inletnumber++;
			myfile << inletnumber << "," << 0.5*(openings[i].innerRadius + openings[i].outerRadius)*dx << "," << openings[i].area*dx*dx << "\n";
		}

	}
	myfile.close();
	pcout << "-> writing Radii file\n";
}


void createLinkFile(
		MultiScalarField3D<VoxelFlagType>& voxelMatrix,
		MultiTensorField3D<plint,159>& linkField,
		std::string fname)
{
	plint X = voxelMatrix.getBoundingBox().x1-voxelMatrix.getBoundingBox().x0;
	plint Y = voxelMatrix.getBoundingBox().y1-voxelMatrix.getBoundingBox().y0;
	plint Z = voxelMatrix.getBoundingBox().z1-voxelMatrix.getBoundingBox().z0;

	MultiScalarField3D<VoxelFlagType> flagMatrix((MultiBlock3D&) voxelMatrix);

	setToConstant(flagMatrix,
			flagMatrix.getBoundingBox(), static_cast<VoxelFlagType>(0));
	setToConstant(flagMatrix, voxelMatrix, voxelFlag::inside,
			flagMatrix.getBoundingBox(), static_cast<VoxelFlagType>(1));
	if (linkFileUseInnerBorder) {
		setToConstant(flagMatrix, voxelMatrix,
				voxelFlag::innerBorder, flagMatrix.getBoundingBox(), static_cast<VoxelFlagType>(1));}
	if (linkFileUseOuterBorder) {
		setToConstant(flagMatrix, voxelMatrix,
				voxelFlag::outerBorder, flagMatrix.getBoundingBox(), static_cast<VoxelFlagType>(1));}

	double numCellsAll = computeDoubleSum(flagMatrix);
	pcout << "--> total number of fluid cells: " << numCellsAll << std::endl;

	{
		// This bounding box is our overall link checking
		Box3D bbox = voxelMatrix.getBoundingBox();

		// Let us make the hash container of the triangles
		// This part is copied from the voxelizer
		MultiContainerBlock3D hashContainer(voxelMatrix);

    	std::vector<MultiBlock3D*> container_arg;
    		container_arg.push_back(&hashContainer);

		// Generate the hashes 
    	applyProcessingFunctional (
            new CreateTriangleHash<T>(*arterySurfaceMesh),
            hashContainer.getBoundingBox(), container_arg );	

		pcout << "--> Applying Processing Functional to compute links...\n";
		double t1 = MPI_Wtime();

		// Like creating the Hash, we need to pack up the arguments into a
		// vector of MultiBlock3D*-s
		std::vector<MultiBlock3D*> voxelToLinkArgs;
    	voxelToLinkArgs.push_back(&flagMatrix);
		voxelToLinkArgs.push_back(&linkField);
		voxelToLinkArgs.push_back(&hashContainer);

		// Thisis where we generate the links
		applyProcessingFunctional(new VoxelToLink3D<VoxelFlagType,plint>(bbox), bbox, voxelToLinkArgs);
		double t2 = MPI_Wtime();
		pcout << "--> done\n";
		pcout << "--> VoxelToLink3D data processor took: " << t2 - t1 << " sec.\n";
	}

}

// this is the function that prepares and performs the voxelization
void run(bool endearly)
{
	// get the initial bounding box before palabos centres origin to 0
	Cuboid<T> cuboid = (*arteryTriangleSet).getBoundingCuboid(); 
	plint resolution = (plint)(std::round((cuboid.x1()-cuboid.x0())/dxREL));

	// Rescale ArteryTriangleSet here
	(*arteryTriangleSet).scale(dxREL/(((cuboid.x1()-cuboid.x0())/(T)(resolution))));
	pcout << "cuboid length for dX to fit: " << std::setprecision(10) << (cuboid.x1()-cuboid.x0()) << std::endl;
	pcout << "Resolution calculated for dX to fit: " << std::setprecision(10) << resolution << " and " << (T)(resolution) << std::endl;
	pcout << "rescale factor to get dX to fit: " << std::setprecision(10) << dxREL/(((cuboid.x1()-cuboid.x0())/(T)(resolution))) << std::endl;

	// Create new DEFMesh with rescaled arteryTriangleSet
	DEFscaledMesh<T>* arteryDefMesh = 
		new DEFscaledMesh<T>(*arteryTriangleSet, resolution, referenceDirection, margin, extraLayer);

	pcout << "-> dx: " << std::setprecision(10) << arteryDefMesh->getDx() << " in units of .stl" << std::endl;
	
    double dx = dxABS;
	double shift_x = -cuboid.x0()/dxREL+3.0, shift_y = -cuboid.y0()/dxREL+3.0, shift_z = -cuboid.z0()/dxREL+3.0;

	TriangleBoundary3D<T>* arteryBoundary=new TriangleBoundary3D<T>(*arteryDefMesh);
	(*arteryBoundary).getMesh().inflate();

	// translate the mesh
	Array<T,2> xRange, yRange, zRange;
	(*arteryBoundary).getMesh().computeBoundingBox(xRange, yRange, zRange);
	(*arteryBoundary).getMesh().translate(Array<T,3>(-xRange[0]+3.0,-yRange[0]+3.0,-zRange[0]+3.0));
	
	arterySurfaceMesh = new TriangularSurfaceMesh<T>((*arteryBoundary).getMesh());

	if( endearly )
	{
		pcout << "-> dumping ArteryTriangleSet.stl" << std::endl;
		TriangleSet<T> ArteryTriangleSet = arterySurfaceMesh->toTriangleSet(DBL);
		ArteryTriangleSet.writeBinarySTL("ArteryTriangleSet.stl");
	}

	pluint num_openings = (plint)(*arteryBoundary).getInletOutlet().size();
	plint sortDirection = 0;
	openings.reserve(num_openings);

	if (num_openings != (pluint)analysisPoints.size()) {
		pcout << "num openings " << num_openings << std::endl;
		pcout << "analysisPoints " << (pluint)analysisPoints.size() << std::endl;
		throw PlbIOException("-> too few analysis points");
	}

	pcout << "-> locating openings and identifying lid triangles" << std::endl;
	for(int i=0; i < analysisPoints.size(); i++) {
		pcout << "Analysis Points[" << i << "]=( " << (analysisPoints[i])[0] 
			<<	" , " << (analysisPoints[i])[1] << " , " << (analysisPoints[i])[2] << " ) \n";
	}
	
	bool amIMainRank = global::mpi().isMainProcessor();
	
	std::ofstream ioletfile;

	if( amIMainRank ) {
		ioletfile.open("ioletpositions.txt");
		ioletfile << "DX: " << dx << std::endl;
		ioletfile << "SHIFTS: " << shift_x << " " << shift_y << " " << shift_z << std::endl;
	}
	
	// identify lid triangles at this opening; store in unordered_map for later
	std::vector<Lid> lidList = (*arteryBoundary).getInletOutlet(sortDirection);

	for (pluint i = 0; i < num_openings; ++i) {
		// set distanceToIOlet to some large value
		T distanceToIOlet = std::numeric_limits<T>::max();

		// compute center and normal of opening
		openings[i].center = computeBaryCenter (
				(*arteryBoundary).getMesh(),
				lidList[i] );
		openings[i].normal = computeNormal (
				(*arteryBoundary).getMesh(),
				lidList[i] );

		// compute radii values
		openings[i].innerRadius = computeInnerRadius (
				(*arteryBoundary).getMesh(),
				lidList[i] );

		openings[i].outerRadius = computeOuterRadius (
				(*arteryBoundary).getMesh(),
				lidList[i] );

		//JM Areas to go here.
		openings[i].area = computeArea (
				(*arteryBoundary).getMesh(),
				lidList[i] );

		// identify opening as inlet or outlet
		for (pluint j = 0; j < num_openings; ++j) {
			if (norm(analysisPoints[j]-openings[i].center) < distanceToIOlet) {
				distanceToIOlet = norm(analysisPoints[j]-openings[i].center);
				if (j < numIlets) openings[i].type = 2;
				else			  openings[i].type = 3;
				openings[i].index = j;
			}
		}

		// output to screen
		pcout << "-> " << i <<" :  opening at  (" << openings[i].center[0] << "," << openings[i].center[1] << ","
			<< openings[i].center[2] << ")" << std::endl;
		pcout << "-----> normal: (" << openings[i].normal[0] << "," << openings[i].normal[1] << ","
			<< openings[i].normal[2] << ")" << std::endl;
		pcout << "------> index: (" << openings[i].index << ")" << std::endl;
		pcout << "-------> type: (" << openings[i].type << ")" << std::endl;

		pluint frstTriangle = lidList[i].firstTriangle;
		pluint numTriangles = lidList[i].numTriangles;
		for (pluint iTriangle = frstTriangle; iTriangle < frstTriangle+numTriangles; ++iTriangle) {
			lidTriangles[iTriangle] = i;
		//	pcout << "Marking " << iTriangle << " as a LID for opening " << i << "\n";
		}
		
		if( amIMainRank ) {
				ioletfile << openings[i].center[0] << " " << openings[i].center[1] << " " << openings[i].center[2] << " " << 0.5*(openings[i].innerRadius + openings[i].outerRadius) << std::endl;
				//ioletfile << openings[i].center[0] << " " << openings[i].center[1] << " " << openings[i].center[2] << std::endl;
		}
			// if (openings[i].innerRadius < 1.0)
		if (0.5*(openings[i].innerRadius + openings[i].outerRadius) < 1.0)
		{
			if( amIMainRank ) ioletfile.close();
			throw PlbIOException("Average Radius smaller than dx at above location (last point in ioletpostions.txt)");
		}
	}

	if ( amIMainRank ) ioletfile.close();

	if(endearly == true) {	
		pcout << "Ending early. Positions were output to ioletpositions.txt" << std::endl;
		// exit(0); //JM Robin's original
		return;
	}

	memu();

	// VOXELIZATION PROCESS
	pcout << "-> voxelizing" << std::endl;
	const VoxelFlagType arteryFlowType = voxelFlag::inside;
 	double start = MPI_Wtime();
	VoxelizedDomain3D<T>* arteryVoxelizedDomain = new VoxelizedDomain3D<T>(
			(*arteryBoundary), arteryFlowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);
	double end = MPI_Wtime();
	memu();
	pcout << "--> Voxelization took " << end -start << " sec.\n";

	if (linkFileUseInnerBorder) {
		if (linkFileUseOuterBorder) {
			pcout << "--> links for inner and outer border will be produced" << std::endl;
		} else {
			pcout << "--> links for inner border will be produced, but not for outer border" << std::endl;
		}
	} else {
		linkFileUseOuterBorder = false;
		pcout << "--> links will not be produced for either inner or outer border" << std::endl;
	}

    memu();
	pcout << "-> Allocating Link field based on voxel matrix\n";
    MultiTensorField3D<plint,159> linkField(arteryVoxelizedDomain->getVoxelMatrix());
	memu();
	pcout << "-> Creating links\n"; 

	start = MPI_Wtime();
	createLinkFile(arteryVoxelizedDomain->getVoxelMatrix(), linkField, "fluidAndLinks.dat");
	end = MPI_Wtime();
	pcout << "-> Links creation took " << end - start << " sec\n";
	memu();
	pcout << "-> Deleting Voxelized Domain\n"; 

	delete arteryVoxelizedDomain;
	delete arteryBoundary;
	delete arteryDefMesh;

	memu();
	std::string linkFileName = "fluidsAndLinks.dat";
	pcout << "-> Outputing Links file to " << linkFileName << "\n";
	linkField.setDataTransfer(new SparseTensorFieldDataTransfer3D<plint,159>());
	start = MPI_Wtime();
	parallelIO::save(linkField, linkFileName.c_str(), false);
	end = MPI_Wtime();
	pcout << "-> ParallelIO took " << end - start << " sec.\n";

	// write the heme input xml file with the correct inlet and outlet positions (in physical units)
	if( amIMainRank ) writeXML(num_openings, dx, shift_x, shift_y, shift_z);
	pcout << "-> done" << std::endl;
	return;
}

int main(int argc, char* argv[])
{
	// end early if two arguments are given instead of one
	bool endearly = (argc == 3);
	
	// Need to init PLB before using pcout otherwise
	// every MPI will write to stdout 
	plbInit(&argc, &argv);

	if(endearly) {
		pcout << "ENDEARLY option: ending early (after iolets positions dump)" << std::endl;
	}

	global::directories().setOutputDir("./");
	global::IOpolicy().activateParallelIO(true);

	std::string generalParamXmlFileName;
	try {
		global::argv(1).read(generalParamXmlFileName);
	}
	catch (PlbIOException& exception) {
		pcout << "-> wrong parameters; the syntax is: "
			<< (std::string) global::argv(0) << " params.xml" << std::endl;
		return -1;
	}

	// read the parameter xml input file
	try {
		XMLreader generalXML(generalParamXmlFileName);
		readGeneralParameters(generalXML);
	}
	catch (PlbIOException& exception) {
		pcout << "-> error in input file " << generalParamXmlFileName
			<< ": " << exception.what() << std::endl;
		return -1;
	}

	try {
		run(endearly);
	}
	catch (PlbIOException& exception) {
		pcout << exception.what() << std::endl;
		return -1;
	}

	return 0;
}
