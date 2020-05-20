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

#ifndef OBSOLETE_FORMAT_WRAPPER_H
#define OBSOLETE_FORMAT_WRAPPER_H

#include "geometry/rawTriangleMesh.h"
#include "offLattice/triangleBoundary3D.h"
#include <string>
#include <map>

namespace plb {


template<typename T>
RawTriangleMesh<T> triangleSetToRawTriangleMesh(TriangleSet<T> const& triangleSet, T eps=getEpsilon<T>(DBL));

template <typename T>
RawConnectedTriangleMesh<T> triangleSetToConnectedTriangleMesh(TriangleSet<T> const& triangleSet, T eps=getEpsilon<T>(DBL));

template <typename T>
TriangleSet<T> rawTriangleMeshToTriangleSet(RawTriangleMesh<T> const& triangleMesh, Precision precision = DBL );

template <typename T>
TriangleSet<T> rawTriangleMeshToTriangleSet(RawTriangleMesh<T> const& triangleMesh, T eps );

template <typename T>
RawConnectedTriangleMesh<T> def_to_ConnectedMesh(TriangleBoundary3D<T>& boundary, std::string partNaming = "Part");

}  // namespace plb

#endif  // OBSOLETE_FORMAT_WRAPPER_H

