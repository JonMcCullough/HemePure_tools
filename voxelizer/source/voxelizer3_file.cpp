#include "palabos3D.h"
#include "palabos3D.hh"

#include <cmath>
#include <cstdio>
#include <limits>
#include <unordered_map>

#include <iostream>
#include <fstream>

#include "utils.h"

// ========== DO NOT TOUCH (START) ========== //

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

const plint blockSize = 0;				// Zero means: no sparse representation.
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

// defining a hash map to store triangles near site
// ------------------------------------------------
struct funcHash {
	size_t operator()(const Dot3D &k) const{
		size_t h1 = std::hash<T>()(k.x);
		size_t h2 = std::hash<T>()(k.y);
		size_t h3 = std::hash<T>()(k.z);
		return (h1^(h2 << 1))^(h3 << 2);
	}
};

struct funcEqul {
	bool operator()(const Dot3D& lhs, const Dot3D& rhs) const{
		return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
	}
};

//typedef std::unordered_map<Dot3D, std::vector<plint>, funcHash, funcEqul> Dot3DMap;
//Dot3DMap triangles;
std::vector<plint> triangles;
// ---

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

template<typename T1, typename T2>
class VoxelToLink3D : public BoxProcessingFunctional3D_ST<T1, T2, 159> {
	public:
		VoxelToLink3D(Box3D bbox);
		virtual void process(Box3D domain, ScalarField3D<T1>& flags, TensorField3D<T2,159>& links);
		virtual VoxelToLink3D<T1,T2>* clone() const {
			return new VoxelToLink3D(*this);
		}
		virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
			modified[0] = modif::nothing;
			modified[1] = modif::staticVariables;
		}
	private:
		T2 globalid(Dot3D const& absPos, plint iX, plint iY, plint iZ) const;
		void getTriangles(
				Array<T2,2> const& xRange,
				Array<T2,2> const& yRange,
				Array<T2,2> const& zRange,
				std::vector<plint>& foundTriangles) const;
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

template<typename T1, typename T2>
T2 VoxelToLink3D<T1,T2>::globalid(
		Dot3D const& absPos, plint iX, plint iY, plint iZ) const {
	return iZ+absPos.z+nZ*(iY+absPos.y+nY*(iX+absPos.x));
}

//template<typename T1, typename T2>
//void VoxelToLink3D<T1,T2>::getTriangles(
//		Array<T2,2> const& xRange,
//		Array<T2,2> const& yRange,
//		Array<T2,2> const& zRange,
//		std::vector<plint>& foundTriangles) const
//{
//	// fit onto the grid by making it bigger
//	// to be sure the triangle is never missed through round-off errors
//	Box3D discreteRange(
//			(plint)xRange[0], (plint)xRange[1]+1,
//			(plint)yRange[0], (plint)yRange[1]+1,
//			(plint)zRange[0], (plint)zRange[1]+1);
//
//	foundTriangles.clear();
//	for (plint iX = discreteRange.x0; iX <= discreteRange.x1; ++iX)
//		for (plint iY = discreteRange.y0; iY <= discreteRange.y1; ++iY)
//			for (plint iZ = discreteRange.z0; iZ <= discreteRange.z1; ++iZ) {
//				if (triangles.find(Dot3D(iX,iY,iZ)) != triangles.end()) {
//					std::vector<plint> const& newTriangles = triangles.at(Dot3D(iX,iY,iZ));
//					foundTriangles.insert(foundTriangles.end(),
//							newTriangles.begin(), newTriangles.end());
//				}
//			}
//
//	std::sort(foundTriangles.begin(), foundTriangles.end());
//	foundTriangles.erase(unique(foundTriangles.begin(), foundTriangles.end()), foundTriangles.end());
//}

template<typename T1, typename T2>
void VoxelToLink3D<T1,T2>::process(
		Box3D domain, ScalarField3D<T1>& flags, TensorField3D<T2,159>& links)
{
	Dot3D absPos = flags.getLocation();
	Dot3D ofs = computeRelativeDisplacement(flags, links);
	//pcout << "--> processing links and finding intersections with triangles on " << getpid() << std::endl;

	// variable declarations
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
  //std::vector<T2> possibleTriangles;

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
	for (plint iX = domain.x0; iX <= domain.x1; ++iX)
		for (plint iY = domain.y0; iY <= domain.y1; ++iY)
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

										Array<T2,2> xRange(
												std::min(point1[0], point2[0]),
												std::max(point1[0], point2[0]));
										Array<T2,2> yRange(
												std::min(point1[1], point2[1]),
												std::max(point1[1], point2[1]));
										Array<T2,2> zRange(
												std::min(point1[2], point2[2]),
												std::max(point1[2], point2[2]));

										//getTriangles(xRange, yRange, zRange, possibleTriangles);

										std::vector<T> crossings;
										for (pluint iPossible = 0; iPossible < triangles.size()/*possibleTriangles.size()*/; ++iPossible) {
											plint iTriangle = triangles/*possibleTriangles*/[iPossible];
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
	pcout << "-| memory use: " << vmrss_kb << " kB" << "\n";
}

void writeXML(const pluint num_openings, double dx, double shift_x, double shift_y, double shift_z)
{
	std::ofstream myfile;
	myfile.open ("iolets_block_inputxml.txt");

	myfile << "  <inlets>\n";

	int inletnumber = 0;
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

	//myfile << "outlets\n";
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

// assigns all triangles beforehand, but with heavy memory usage
//void assignTriangles(
//		TriangularSurfaceMesh<T> const& mesh,
//		//std::vector<Dot3D>& assignedPositions,
//		Dot3DMap& triangles)
//{
//	//assignedPositions.clear();
//	for (plint iTriangle = 0; iTriangle < mesh.getNumTriangles(); ++iTriangle)
//	{
//		Array<T,3> const& vertex0 = mesh.getVertex(iTriangle,0);
//		Array<T,3> const& vertex1 = mesh.getVertex(iTriangle,1);
//		Array<T,3> const& vertex2 = mesh.getVertex(iTriangle,2);
//
//		Array<T,2> xRange(
//				std::min(vertex0[0], std::min(vertex1[0], vertex2[0])),
//				std::max(vertex0[0], std::max(vertex1[0], vertex2[0])));
//		Array<T,2> yRange(
//				std::min(vertex0[1], std::min(vertex1[1], vertex2[1])),
//				std::max(vertex0[1], std::max(vertex1[1], vertex2[1])));
//		Array<T,2> zRange(
//				std::min(vertex0[2], std::min(vertex1[2], vertex2[2])),
//				std::max(vertex0[2], std::max(vertex1[2], vertex2[2])));
//
//		// fit onto the grid by making it bigger
//		// to be sure the triangle is never missed through round-off errors
//		Box3D discreteRange(
//				(plint)xRange[0], (plint)xRange[1]+1,
//				(plint)yRange[0], (plint)yRange[1]+1,
//				(plint)zRange[0], (plint)zRange[1]+1);
//
//		for (plint iX = discreteRange.x0; iX <= discreteRange.x1; ++iX)
//			for (plint iY = discreteRange.y0; iY <= discreteRange.y1; ++iY)
//				for (plint iZ = discreteRange.z0; iZ <= discreteRange.z1; ++iZ) {
//					//if (triangles.find(Dot3D(iX,iY,iZ)) == triangles.end()) {
//					//	assignedPositions.push_back(Dot3D(iX,iY,iZ));
//					//}
//					triangles[Dot3D(iX,iY,iZ)].push_back(iTriangle);
//				}
//	}
//}

//void assignTriangles(
//		TriangularSurfaceMesh<T> const& mesh,
//		Box3D const& subdomain,
//		Dot3DMap& triangles)
void assignTriangles(
		Box3D subdomain)
{
	triangles.clear(); pluint hit = 0;
	for (plint iTriangle = 0; iTriangle < arterySurfaceMesh->getNumTriangles(); ++iTriangle)
	{
		Array<T,3> const& vertex0 = arterySurfaceMesh->getVertex(iTriangle,0);
		Array<T,3> const& vertex1 = arterySurfaceMesh->getVertex(iTriangle,1);
		Array<T,3> const& vertex2 = arterySurfaceMesh->getVertex(iTriangle,2);

		Array<T,2> xRange(
				std::min(vertex0[0], std::min(vertex1[0], vertex2[0])),
				std::max(vertex0[0], std::max(vertex1[0], vertex2[0])));
		Array<T,2> yRange(
				std::min(vertex0[1], std::min(vertex1[1], vertex2[1])),
				std::max(vertex0[1], std::max(vertex1[1], vertex2[1])));
		Array<T,2> zRange(
				std::min(vertex0[2], std::min(vertex1[2], vertex2[2])),
				std::max(vertex0[2], std::max(vertex1[2], vertex2[2])));

		// fit onto the grid by making it bigger
		// to be sure the triangle is never missed through round-off errors
		Box3D discreteRange(
				(plint)xRange[0], (plint)xRange[1]+1,
				(plint)yRange[0], (plint)yRange[1]+1,
				(plint)zRange[0], (plint)zRange[1]+1);

		subdomain.enlarge(1);

		Box3D inters;
		if (intersect(discreteRange, subdomain, inters)) {
			//for (plint iX = subdomain.x0; iX <= subdomain.x1; ++iX)
			//	for (plint iY = subdomain.y0; iY <= subdomain.y1; ++iY)
			//		for (plint iZ = subdomain.z0; iZ <= subdomain.z1; ++iZ) {
			//			triangles[Dot3D(iX,iY,iZ)].push_back(iTriangle);
			//			hit++;
			//		}
			triangles.push_back(iTriangle);
		} //pcout << "--> number of hits: " << hit << ", with triangles.size() = " << triangles.size() << std::endl;
	}
}

void createLinkFile(
		MultiScalarField3D<int>& voxelMatrix,
		MultiContainerBlock3D& triangleHash,
		std::string fname)
{
	plint X = voxelMatrix.getBoundingBox().x1-voxelMatrix.getBoundingBox().x0;
	plint Y = voxelMatrix.getBoundingBox().y1-voxelMatrix.getBoundingBox().y0;
	plint Z = voxelMatrix.getBoundingBox().z1-voxelMatrix.getBoundingBox().z0;

	MultiScalarField3D<int> flagMatrix((MultiBlock3D&) voxelMatrix);

	setToConstant(flagMatrix,
			flagMatrix.getBoundingBox(), 0);
	setToConstant(flagMatrix, voxelMatrix, voxelFlag::inside,
			flagMatrix.getBoundingBox(), 1);
	if (linkFileUseInnerBorder) {
		setToConstant(flagMatrix, voxelMatrix,
				voxelFlag::innerBorder, flagMatrix.getBoundingBox(), 1);}
	if (linkFileUseOuterBorder) {
		setToConstant(flagMatrix, voxelMatrix,
				voxelFlag::outerBorder, flagMatrix.getBoundingBox(), 1);}

	pluint numCellsAll = computeSum(flagMatrix);
	pcout << "--> total number of fluid cells: " << numCellsAll << std::endl;

	//memu();
	////std::vector<Dot3D> assignedPositions;
	//assignTriangles(*arterySurfaceMesh, /*assignedPositions, */triangles);
	//memu();

	pluint count = 0, progress = 0, numCells = 0;
	pluint maxbb = 100;
	T percent;
	for (plint x0 = voxelMatrix.getBoundingBox().x0; x0 < X; x0 += (maxbb+1))
		for (plint y0 = voxelMatrix.getBoundingBox().y0; y0 < Y; y0 += (maxbb+1))
			for (plint z0 = voxelMatrix.getBoundingBox().z0; z0 < Z; z0 += (maxbb+1)) {

				plint x1 = x0+maxbb; plint y1 = y0+maxbb; plint z1 = z0+maxbb;

				if (x1 >= X) x1 = X;
				if (y1 >= Y) y1 = Y;
				if (z1 >= Z) z1 = Z;

				Box3D bbox(x0,x1,y0,y1,z0,z1);
				numCells = computeSum(flagMatrix, bbox); progress += numCells;

				if (numCells != 0) {
					percent = (T)progress/numCellsAll*100;
					pcout << "--> " <<
						std::setiosflags(std::ios::scientific) << std::setw(5) <<
						percent << "% number of fluid cells: " << numCells << std::endl;

				  //assignTriangles(*arterySurfaceMesh, bbox, triangles);
					assignTriangles(bbox);

					memu();
					MultiTensorField3D<plint,159> linkField(voxelMatrix, bbox);
					memu();

					char buffer[32]; count++;
					snprintf(buffer, sizeof(char)*32, "fluidAndLinks_%03li.dat", count);

					pcout << "--> " <<
						std::setiosflags(std::ios::scientific) << std::setw(5) <<
						percent << "% processing links and finding intersections with triangles" << std::endl;
					applyProcessingFunctional(new VoxelToLink3D<int,plint>(bbox), bbox, flagMatrix, linkField);
					bool dynamicContent = false;
					pcout << "--> " <<
						std::setiosflags(std::ios::scientific) << std::setw(5) <<
						percent << "% data transfer" << std::endl;
					linkField.setDataTransfer(new SparseTensorFieldDataTransfer3D<plint,159>());
					pcout << "--> " <<
						std::setiosflags(std::ios::scientific) << std::setw(5) <<
						percent << "% parallel save" << std::endl;
					parallelIO::save(linkField, buffer, dynamicContent);
				}
			}
}

// this is the function that prepares and performs the voxelization
void run(bool endearly)
{
	Cuboid<T> cuboid = (*arteryTriangleSet).getBoundingCuboid(); // get the initial bounding box before palabos centres origin to 0
	
	plint resolution = (plint)(std::round((cuboid.x1()-cuboid.x0())/dxREL));

	DEFscaledMesh<T>* arteryDefMesh =
		new DEFscaledMesh<T>(*arteryTriangleSet, resolution, referenceDirection, margin, extraLayer);
	pcout << "-> dx: " << std::setprecision(10) << arteryDefMesh->getDx() << " in units of .stl" << std::endl;

	(*arteryTriangleSet).scale(dxREL/(((cuboid.x1()-cuboid.x0())/(T)(resolution))));
	pcout << "cuboid length for dX to fit: " << std::setprecision(10) << (cuboid.x1()-cuboid.x0()) << std::endl;
	pcout << "Resolution calculated for dX to fit: " << std::setprecision(10) << resolution << "and " << (T)(resolution) << std::endl;
	pcout << "rescale factor to get dX to fit: " << std::setprecision(10) << dxREL/(((cuboid.x1()-cuboid.x0())/(T)(resolution))) << std::endl;

	

	arteryDefMesh =
		new DEFscaledMesh<T>(*arteryTriangleSet, resolution, referenceDirection, margin, extraLayer);
	TriangleBoundary3D<T> arteryBoundary(*arteryDefMesh);
	pcout << "-> dx: " << std::setprecision(10) << arteryDefMesh->getDx() << " in units of .stl" << std::endl;
	delete arteryDefMesh;
	arteryBoundary.getMesh().inflate();

	// translate the mesh
	Array<T,2> xRange, yRange, zRange;
	arteryBoundary.getMesh().computeBoundingBox(xRange, yRange, zRange);
	arteryBoundary.getMesh().translate(Array<T,3>(-xRange[0]+3.0,-yRange[0]+3.0,-zRange[0]+3.0));

	// For reverse transform later
	//double dx = arteryDefMesh->getDx();
	double dx = dxABS;
	double shift_x = -cuboid.x0()/dxREL+3.0, shift_y = -cuboid.y0()/dxREL+3.0, shift_z = -cuboid.z0()/dxREL+3.0;

	arterySurfaceMesh = new TriangularSurfaceMesh<T>(arteryBoundary.getMesh());

	pcout << "-> dumping ArteryTriangleSet.stl" << std::endl;
	TriangleSet<T> ArteryTriangleSet = arterySurfaceMesh->toTriangleSet(DBL);
	ArteryTriangleSet.writeBinarySTL("ArteryTriangleSet.stl");

	pluint num_openings = (plint)arteryBoundary.getInletOutlet().size();
	plint sortDirection = 0;
	openings.reserve(num_openings);

	if (num_openings != (pluint)analysisPoints.size()) {
		throw PlbIOException("-> too few analysis points");
	}

	pcout << "-> locating openings and identifying lid triangles" << std::endl;

	std::ofstream ioletfile;
	ioletfile.open("ioletpositions.txt");
	ioletfile << "DX: " << dx << std::endl;
	ioletfile << "SHIFTS: " << shift_x << " " << shift_y << " " << shift_z << std::endl;
	for (pluint i = 0; i < num_openings; ++i) {
		// set distanceToIOlet to some large value
		T distanceToIOlet = std::numeric_limits<T>::max();

		// compute center and normal of opening
		openings[i].center = computeBaryCenter (
				arteryBoundary.getMesh(),
				arteryBoundary.getInletOutlet(sortDirection)[i] );
		openings[i].normal = computeNormal (
				arteryBoundary.getMesh(),
				arteryBoundary.getInletOutlet(sortDirection)[i] );

		// compute radii values
		openings[i].innerRadius = computeInnerRadius (
				arteryBoundary.getMesh(),
				arteryBoundary.getInletOutlet(sortDirection)[i] );
		openings[i].outerRadius = computeOuterRadius (
				arteryBoundary.getMesh(),
				arteryBoundary.getInletOutlet(sortDirection)[i] );
		//JM Areas to go here.
		openings[i].area = computeArea (
				arteryBoundary.getMesh(),
				arteryBoundary.getInletOutlet(sortDirection)[i] );

		// identify opening as inlet or outlet
		for (pluint j = 0; j < num_openings; ++j) {
			if (norm(analysisPoints[j]-openings[i].center) < distanceToIOlet) {
				distanceToIOlet = norm(analysisPoints[j]-openings[i].center);
				if (j < numIlets) openings[i].type = 2;
				else			  openings[i].type = 3; openings[i].index = j;
			}
		}

		// identify lid triangles at this opening; store in unordered_map for later
		pluint frstTriangle = arteryBoundary.getInletOutlet(sortDirection)[i].firstTriangle;
		pluint numTriangles = arteryBoundary.getInletOutlet(sortDirection)[i].numTriangles;
		for (pluint iTriangle = frstTriangle; iTriangle < frstTriangle+numTriangles; ++iTriangle)
			lidTriangles[iTriangle] = i;

		// output to screen
		pcout << "-> opening at  (" << openings[i].center[0] << "," << openings[i].center[1] << ","
			<< openings[i].center[2] << ")" << std::endl;
		pcout << "-----> normal: (" << openings[i].normal[0] << "," << openings[i].normal[1] << ","
			<< openings[i].normal[2] << ")" << std::endl;
		pcout << "------> index: (" << openings[i].index << ")" << std::endl;
		pcout << "-------> type: (" << openings[i].type << ")" << std::endl;

		ioletfile << openings[i].center[0] << " " << openings[i].center[1] << " " << openings[i].center[2] << std::endl;
	}
	ioletfile.close();

	if(endearly == true) {
		pcout << "Ending early. Positions were output to ioletpositions.txt" << std::endl;
		exit(0);
	}

	memu();
	// VOXELIZATION PROCESS
	pcout << "-> voxelizing" << std::endl;
	const int arteryFlowType = voxelFlag::inside;
	VoxelizedDomain3D<T> arteryVoxelizedDomain(
			arteryBoundary, arteryFlowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);
	memu();

	pcout << "-> saving links" << std::endl;
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
	createLinkFile(arteryVoxelizedDomain.getVoxelMatrix(), arteryVoxelizedDomain.getTriangleHash(), "fluidAndLinks.dat");

	// write the heme input xml file with the correct inlet and outlet positions (in physical units)
	writeXML(num_openings, dx, shift_x, shift_y, shift_z);
	pcout << "-> done" << std::endl;
	return;
}

int main(int argc, char* argv[])
{
	// end early if two arguments are given instead of one
	bool endearly = false;
	pcout << argc << std::endl;
	if(argc == 3) {
		pcout << "2 arguments given - ending early (after iolets positions dump)";
		endearly = true;
	}

	plbInit(&argc, &argv);
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
//pcout << "---> DEBUG" << std::endl;
//{
//	int i = 0;
//	char hostname[256];
//	gethostname(hostname, sizeof(hostname));
//	printf("PID %d on %s ready for attach\n", getpid(), hostname);
//	fflush(stdout);
//	while (0 == i)
//		sleep(5);
//}
