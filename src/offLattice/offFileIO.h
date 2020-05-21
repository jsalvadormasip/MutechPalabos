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

#ifndef OFF_FILE_IO_H
#define OFF_FILE_IO_H

/*
 * OFF I/O stuff. Everything is independent of whether it's going to be used
 * on a RawTriangleMesh, or a RawConnectedTriangleMesh, or a DEFtriangleMesh.
 */

#include "core/array.h"
#include "core/globalDefs.h"
#include "offLattice/triangleMesh.h"

#include <vector>

namespace plb {

template<typename T>
class OFFreader {
public:
    OFFreader(std::string fname);
    std::vector<Array<T,3> > const& getVertices() const { return vertices; }
    std::vector<std::vector<plint> > const& getFacets() const { return facets; }
private:
    void readOFF(std::string fname);
    // Caution: in the readAsciiOFF method we keep all information that resides
    //          in the OFF file. In the readAsciiOFF method of the TriangleSet class
    //          we neglect all triangles that have one or more edges with length
    //          equal to 0. We do not do this here, because from the OFFreader we
    //          can directly create a connected mesh and we do not want to break
    //          vertex and triangle numbering.
    bool readAsciiOFF(FILE* fp);
    void skipLines(plint nLines, FILE* fp) const;
    int readAhead(FILE *fp, char commentCharacter) const;
private:
    std::vector<Array<T,3> > vertices;
    std::vector<std::vector<plint> > facets;
};

template<typename T>
void writeAsciiOFF(TriangleMesh<T>& mesh, std::string fname, T eps=getEpsilon<T>(DBL), int numDecimalDigits=10);

}  // namespace plb

#endif  // OFF_FILE_IO_H

