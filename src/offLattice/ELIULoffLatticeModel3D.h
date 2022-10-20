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

#ifndef ELIUL_OFF_LATTICE_MODEL_3D_H
#define ELIUL_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/offLatticeModel3D.h"

namespace plb {

/**
 * See (Marson et al., 2021)
 *
 * (Marson et al., 2021) Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., Latt, J., 2021.
 * Enhanced single-node lattice Boltzmann boundary condition for fluid flows. Phys. Rev. E 103, 053308.
 * https://doi.org/10.1103/PhysRevE.103.053308
 *
 * @tparam T
 * @tparam Descriptor
 */
template <typename T, template <typename U> class Descriptor>
class ELIULlocalModel3D : public OffLatticeModel3D<T, Array<T, 3> > {
public:
    enum class ParametrizationType{no_param=0,k2=2};
    using PT = ParametrizationType;
    enum class NonEquilibriumType{N=1};
    using NEQ = NonEquilibriumType;
    ELIULlocalModel3D(BoundaryShape3D<T,Array<T,3> >* shape_,
                      int flowType_,ParametrizationType param_ = PT::k2,NonEquilibriumType neq_t_ = NEQ::N);
    virtual ELIULlocalModel3D<T, Descriptor> *clone() const;
    virtual plint getNumNeighbors() const;
    virtual bool isExtrapolated() const;
    virtual void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container);
    virtual void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args);
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &boundaryNode,
        std::vector<int> const &solidDirections, std::vector<plint> const &boundaryIds,
        std::vector<bool> const &hasFluidNeighbor, Dot3D const &absoluteOffset,
        Array<T, 3> &localForce, std::vector<AtomicBlock3D *> const &args);
    virtual ContainerBlockData *generateOffLatticeInfo() const;
    virtual Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const;

private:
    std::vector<T> invAB;
    ParametrizationType param;
    NonEquilibriumType neq_t;
private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class ELIULOffLatticeInfo3D : public ContainerBlockData {
    public:
        ELIULOffLatticeInfo3D() : localForce(Array<T, 3>::zero()) { }
        std::vector<Dot3D> const &getBoundaryNodes() const
        {
            return boundaryNodes;
        }
        std::vector<Dot3D> &getBoundaryNodes()
        {
            return boundaryNodes;
        }
        std::vector<std::vector<int> > const &getSolidDirections() const
        {
            return solidDirections;
        }
        std::vector<std::vector<int> > &getSolidDirections()
        {
            return solidDirections;
        }
        std::vector<std::vector<plint> > const &getBoundaryIds() const
        {
            return boundaryIds;
        }
        std::vector<std::vector<plint> > &getBoundaryIds()
        {
            return boundaryIds;
        }
        std::vector<std::vector<bool> > const &getHasFluidNeighbor() const
        {
            return hasFluidNeighbor;
        }
        std::vector<std::vector<bool> > &getHasFluidNeighbor()
        {
            return hasFluidNeighbor;
        }
        Array<T, 3> const &getLocalForce() const
        {
            return localForce;
        }
        Array<T, 3> &getLocalForce()
        {
            return localForce;
        }
        virtual ELIULOffLatticeInfo3D *clone() const
        {
            return new ELIULOffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> boundaryNodes;
        std::vector<std::vector<int> > solidDirections;
        std::vector<std::vector<plint> > boundaryIds;
        std::vector<std::vector<bool> > hasFluidNeighbor;
        Array<T, 3> localForce;
    };
};

}  // namespace plb

#endif  // ELIUL_OFF_LATTICE_MODEL_3D_H

