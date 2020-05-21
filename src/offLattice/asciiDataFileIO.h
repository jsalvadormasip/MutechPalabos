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

#ifndef ASCII_DATA_FILE_IO_H
#define ASCII_DATA_FILE_IO_H

/* ASCII data I/O stuff. Everything is independent of whether it's going to be used
 * on a RawTriangle, or a RawConnectedTriangleMesh, or a DEFtriangleMesh.
 */


#include "core/globalDefs.h"
#include "offLattice/connectedTriangleMesh.h"

namespace plb {

template<typename T>
void writeAsciiData(RawConnectedTriangleMesh<T>& mesh, FileName fname);

}  // namespace plb

#endif  // ASCII_DATA_FILE_IO_H

