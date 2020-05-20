/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ASCII_DATA_FILE_IO_HH
#define ASCII_DATA_FILE_IO_HH

#include <limits>
#include <map>
#include <fstream>
#include "geometry/asciiDataFileIO.h"

namespace plb {

template<typename T>
void writeAsciiData(RawConnectedTriangleMesh<T>& mesh, FileName fname)
{
    typedef typename ConnectedTriangleMesh<T>::PTriangleIterator PTriangleIterator;
    typedef typename ConnectedTriangleMesh<T>::PVertexIterator PVertexIterator;
    typedef typename ConnectedTriangleMesh<T>::PTriangle PTriangle;
    typedef typename ConnectedTriangleMesh<T>::PVertex PVertex;

    if (!global::mpi().isMainProcessor())
        return;

    if (mesh.getNumTriangles() == 0)
        return;

    for (plint iPart=0; iPart<mesh.numParts(); ++iPart) {
        FileName partFileName(fname);
        if (mesh.numParts() > 1) {
            partFileName.setName(createFileName(partFileName.getName()+"_", iPart, 3));
        }
        RawConnectedTriangleMesh<T> part = extractConnectedPart(mesh, iPart);
        std::ofstream ofile(partFileName.get().c_str());
        ofile.precision(10);
        std::scientific(ofile);

        plint vertexIDtagging = part.getVertexTag("UniqueID");
        std::map<plint,plint> toNewVertexID;
        PVertexIterator vertexIt = part.vertexIterator();
        plint newVertexID=0;
        while (!vertexIt->end()) {
            PVertex vertex(vertexIt->next());
            toNewVertexID[vertex->tag(vertexIDtagging)] = newVertexID;
            ++newVertexID;
        }

        ofile << "Surface data created with Palabos\n";
        ofile << "\n";
        ofile << "Number of triangles: " << part.getNumTriangles() << "\n";
        ofile << "Number of vertices: " << part.getNumVertices() << "\n";
        ofile << "\n";
        ofile << "Vertex coordinates:\n";
        vertexIt = part.vertexIterator();
        while (!vertexIt->end()) {
            PVertex vertex = vertexIt->next();
            ofile << (*vertex)[0] << " " << (*vertex)[1] << " " << (*vertex)[2] << "\n";
        }
        ofile << "\n";
        ofile << "Triangles (vertex indices):\n";
        PTriangleIterator triangleIt = part.triangleIterator();
        while (!triangleIt->end()) {
            PTriangle triangle = triangleIt->next();
            plint i0 = triangle->vertex(0)->tag(vertexIDtagging);
            plint i1 = triangle->vertex(1)->tag(vertexIDtagging);
            plint i2 = triangle->vertex(2)->tag(vertexIDtagging);
            ofile << toNewVertexID[i0] << " " << toNewVertexID[i1] << " " << toNewVertexID[i2] << "\n";
        }
        ofile << "\n";
        for (plint property=0; property<mesh.numVertexProperties(); ++property) {
            ofile << "Scalar property: " << mesh.getVertexPropertyName(property) << "\n";
            vertexIt = part.vertexIterator();
            while (!vertexIt->end()) {
                PVertex vertex = vertexIt->next();
                ofile << vertex->property(property) << "\n";
            }
            if (property != mesh.numVertexProperties()-1) {
                ofile << "\n";
            }
        }
    }
}

}  // namespace plb


#endif  // ASCII_DATA_FILE_IO_HH

