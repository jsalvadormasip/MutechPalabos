/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 * 
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at 
 * <https://palabos.unige.ch/>
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

/** \file
 * Interface to NumPy -- header file.
 */
#ifndef BLOCK_NUMPY_INTERFACE_3D_H
#define BLOCK_NUMPY_INTERFACE_3D_H

#include "multiBlock/multiDataField3D.h"

// http://blog.dhananjaynene.com/2009/03/constructor-method-overloading-in-python/

namespace plb {

template<typename T>
class NTensorField2NumPy3D {
public:
    NTensorField2NumPy3D(MultiNTensorField3D<T>& field_);
    NTensorField2NumPy3D(MultiNTensorField3D<T>& field_, Box3D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiNTensorField3D<T>& field;
    Box3D domain;
};

template<typename T>
class NumPy2NTensorField3D {
public:
    NumPy2NTensorField3D(MultiNTensorField3D<T>& field_);
    NumPy2NTensorField3D(MultiNTensorField3D<T>& field_, Box3D const& domain_);
    void execute(T* array, int size);
    int getSize() const;
private:
    MultiNTensorField3D<T>& field;
    Box3D domain;
};

}  // namespace plb

#endif  // NUMPY_INTERFACE_3D_H