#pragma once
#include <cerrno>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define Distance2(x, y)                                                        \
    ((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1])             \
        + (x[2] - y[2]) * (x[2] - y[2]))

namespace plb
{

template <typename T>
void multiProcWriteVTK(
    // mesh in lattice units since it comes from palabos
	RawConnectedTriangleMesh<T>& mesh, FileName fname,
	// Visualization purposes
	pluint nx, pluint ny, pluint nz,
	// Conversion from Lattice to Physical units for visualization
	T dx, T dt, T rho,
	// Legacy reasons
	bool writeVertexNormals = false, std::string vertexNormalsName = ""
    )
{
	if (mesh.getNumTriangles() == 0)
		return;

	typedef typename ConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
	typedef typename ConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
	typedef typename ConnectedTriangleMesh<T>::PTriangle PTriangle;
	typedef typename ConnectedTriangleMesh<T>::PVertex PVertex;
		
	std::ofstream ofile;
	ofile.open(fname.get().c_str(), std::ofstream::out);
	ofile.precision(10);
	std::scientific(ofile);

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Surface mesh created with Palabos\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";
	ofile << "POINTS " << mesh.getNumVertices() << (sizeof(T) == sizeof(double) ? " double" : " float") << "\n";

	plint vertexIDtagging = mesh.getVertexTag("UniqueID");
	PVertexIterator vertexIt = mesh.vertexIterator();
	std::map<plint, plint> toNewVertexID;
	plint newVertexID = 0;
	while (!vertexIt->end())
    {
		PVertex vertex = vertexIt->next();
		// From lattice to physical
		ofile << (*vertex)[0] * dx << " " << (*vertex)[1] * dx << " " << (*vertex)[2] * dx << "\n";
		toNewVertexID[vertex->tag(vertexIDtagging)] = newVertexID;
		++newVertexID;
	}
	ofile << "\n";

	T n = (T)MIN(nx, MIN(ny, nz)) / 3.;
	plint triangles_to_render = 0;
	PTriangleIterator triangleIt = mesh.triangleIterator();
	//bool RenderBody = true;
	while (!triangleIt->end())
    {
		PTriangle triangle = triangleIt->next();
		//plint i0 = triangle->vertex(0)->tag(vertexIDtagging);
		//plint i1 = triangle->vertex(1)->tag(vertexIDtagging);
		//plint i2 = triangle->vertex(2)->tag(vertexIDtagging);
		Array<T, 3> v0 = triangle->vertex(0)->get();
		Array<T, 3> v1 = triangle->vertex(1)->get();
		Array<T, 3> v2 = triangle->vertex(2)->get();
		T d0 = sqrt(Distance2(v0, v1));
		T d1 = sqrt(Distance2(v0, v2));
		T d2 = sqrt(Distance2(v2, v1));
		
        if (d0 < n && d1 < n && d2 < n)
			++triangles_to_render;
		//else
		//  RenderBody = false;
	}

	ofile << "CELLS " << triangles_to_render << " " << 4 * triangles_to_render << "\n";

	PTriangleIterator triangleIt_new = mesh.triangleIterator();
	while (!triangleIt_new->end())
    {
		PTriangle triangle = triangleIt_new->next();
		plint i0 = triangle->vertex(0)->tag(vertexIDtagging);
		plint i1 = triangle->vertex(1)->tag(vertexIDtagging);
		plint i2 = triangle->vertex(2)->tag(vertexIDtagging);
		Array<T, 3> v0 = triangle->vertex(0)->get();
		Array<T, 3> v1 = triangle->vertex(1)->get();
		Array<T, 3> v2 = triangle->vertex(2)->get();
		T d0 = sqrt(Distance2(v0, v1));
		T d1 = sqrt(Distance2(v0, v2));
		T d2 = sqrt(Distance2(v2, v1));
		if (d0 < n && d1 < n && d2 < n)
		{
			//if (RenderBody)
			ofile << "3 " << i0 << " " << i1 << " " << i2 << "\n";
		}
	}
	ofile << "\n";

	ofile << "CELL_TYPES " << triangles_to_render << "\n";

	for (plint i = 0; i < triangles_to_render; ++i)
	    ofile << "5\n";
	
	ofile << "\n";

    // Do not need any additional data linked to the mesh
    /*
	ofile << "POINT_DATA " << mesh.getNumVertices() << "\n";

	// Only for forces conversion from lattice to physical
	T Cf = rho * (dx * dx * dx * dx) / (dt * dt);
	for (plint property = 0; property < mesh.numVertexProperties(); ++property) {
		ofile << "SCALARS " << mesh.getVertexPropertyName(property)
				<< (sizeof(T) == sizeof(double) ? " double" : " float") << " 1\n"
				<< "LOOKUP_TABLE default\n";
		vertexIt = mesh.vertexIterator();
		while (!vertexIt->end()) {
			PVertex vertex = vertexIt->next();
			// from lattice to physical
			ofile << vertex->property(property) * Cf << "\n";
		}
		ofile << "\n";
	}

	if (writeVertexNormals) {
		ofile << "VECTORS " << vertexNormalsName
				<< (sizeof(T) == sizeof(double) ? " double" : " float") << "\n";
		vertexIt = mesh.vertexIterator();
		while (!vertexIt->end()) {
			PVertex vertex = vertexIt->next();
			Array<T, 3> n = vertex->normal();
			ofile << n[0] << " " << n[1] << " " << n[2] << "\n";
		}
		ofile << "\n";
	}
    */

	ofile.close();
}
	
}

/*
template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice, plint iter)
{
  const plint imSize = 600;
  const plint nx = lattice.getNx();
  const plint ny = lattice.getNy();
  const plint nz = lattice.getNz();
  const plint zComponent = 2;

  Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
  ImageWriter<T> imageWriter("leeloo");

  imageWriter.writeScaledGif( createFileName("uz", iter, 6),
                              *computeVelocityComponent (lattice, slice,
zComponent),
                              imSize, imSize );

  imageWriter.writeScaledGif( createFileName("uNorm", iter, 6),
                              *computeVelocityNorm (lattice, slice),
                              imSize, imSize );
  imageWriter.writeScaledGif( createFileName("omega", iter, 6),
                              *computeNorm(*computeVorticity (
                                      *computeVelocity(lattice) ), slice ),
                              imSize, imSize );
}
*/
