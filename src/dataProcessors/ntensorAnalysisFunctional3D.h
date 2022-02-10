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
 * Helper functions for domain initialization -- header file.
 */
#ifndef NTENSOR_ANALYSIS_FUNCTIONAL_3D_H
#define NTENSOR_ANALYSIS_FUNCTIONAL_3D_H

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template <typename T1, typename T2>
class CopyConvertNTensorFunctional3D : public BoxProcessingFunctional3D_NN<T1, T2> {
public:
    virtual void process(Box3D domain, NTensorField3D<T1> &field1, NTensorField3D<T2> &field2);
    virtual CopyConvertNTensorFunctional3D<T1, T2> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ConvertNTensorToScalarFunctional3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ConvertNTensorToScalarFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template <typename T>
class ConvertScalarToNTensorFunctional3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D *> fields);
    virtual ConvertScalarToNTensorFunctional3D<T> *clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT> &modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // NTENSOR_ANALYSIS_FUNCTIONAL_3D_H
