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

#ifndef ELI_OFF_LATTICE_MODEL_3D_H
#define ELI_OFF_LATTICE_MODEL_3D_H

#include "core/globalDefs.h"
#include "offLattice/guoOffLatticeModel3D.h"
#include "offLattice/offLatticeModel3D.h"
#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus > 201703L)
#include <concepts>
#endif

namespace plb {

/**
 * This class implements the ELI infinite class of directional (link-wise) schemes.
 * Some schemes are named, e.g. ELIUL, ELIULT, ELIFL, others are not and can defined by explicitly
 * providing the alpha- and K- parameters to the constructor using the class LIgeneric.
 *
 * @tparam T
 * @tparam Descriptor
 */
template <typename T, template <typename U> class Descriptor>
class ELIModels3D : public OffLatticeModel3D<T, Array<T, 3>> {
public:
    ELIModels3D(BoundaryShape3D<T, Array<T, 3>> *shape_, int flowType_);
    ELIModels3D<T, Descriptor> *clone() const override = 0;
    [[nodiscard]] plint getNumNeighbors() const override;
    [[nodiscard]] bool isExtrapolated() const override;
    void prepareCell(Dot3D const &cellLocation, AtomicContainerBlock3D &container) override;
    void boundaryCompletion(
        AtomicBlock3D &lattice, AtomicContainerBlock3D &container,
        std::vector<AtomicBlock3D *> const &args) override;

    ContainerBlockData *generateOffLatticeInfo() const override;
    Array<T, 3> getLocalForce(AtomicContainerBlock3D &container) const override;
    void selectComputeStat(bool flag)
    {
        computeStat = flag;
    }
    bool computesStat() const
    {
        return computeStat;
    }

private:
    void cellCompletion(
        BlockLattice3D<T, Descriptor> &lattice, Dot3D const &guoNode,
        std::vector<int> const &dryNodeFluidDirections, std::vector<plint> const &dryNodeIds,
        Dot3D const &absoluteOffset, Array<T, 3> &localForce,
        std::vector<AtomicBlock3D *> const &args);

private:
    bool computeStat;

private:
    /// Store the location of wall nodes, as well as the pattern of missing vs. known
    ///   populations.
    class OffLatticeInfo3D : public ContainerBlockData {
    public:
        OffLatticeInfo3D() : localForce(Array<T, 3>::zero()) { }
        [[nodiscard]] std::vector<Dot3D> const &getDryNodes() const
        {
            return dryNodes;
        }
        std::vector<Dot3D> &getDryNodes()
        {
            return dryNodes;
        }
        [[nodiscard]] std::vector<std::vector<int>> const &getDryNodeFluidDirections() const
        {
            return dryNodeFluidDirections;
        }
        std::vector<std::vector<int>> &getDryNodeFluidDirections()
        {
            return dryNodeFluidDirections;
        }
        [[nodiscard]] std::vector<std::vector<plint>> const &getDryNodeIds() const
        {
            return dryNodeIds;
        }
        std::vector<std::vector<plint>> &getDryNodeIds()
        {
            return dryNodeIds;
        }
        Array<T, 3> const &getLocalForce() const
        {
            return localForce;
        }
        Array<T, 3> &getLocalForce()
        {
            return localForce;
        }
        OffLatticeInfo3D *clone() const override
        {
            return new OffLatticeInfo3D(*this);
        }

    private:
        std::vector<Dot3D> dryNodes;
        std::vector<std::vector<int>> dryNodeFluidDirections;
        std::vector<std::vector<plint>> dryNodeIds;
        Array<T, 3> localForce;
    };
    std::tuple<T, Array<T, 3>> getRhoBarJ(
        const BlockLattice3D<T, Descriptor> &lattice, const Dot3D &guoNode,
        const std::vector<AtomicBlock3D *> &args, const Dot3D &fluidDirection,
        const Cell<T, Descriptor> &cellF) const;

    virtual inline std::array<T, 5> eliCoefficients(T q, T tauPlus, T tauMinus) const = 0;
};

template <typename T, template <typename U> class D>
class ELIUL : public ELIModels3D<T, D> {
    using ELIModels3D<T, D>::ELIModels3D;
    ELIUL<T, D> *clone() const override;
    inline std::array<T, 5> eliCoefficients(T q, T tauPlus, T tauMinus) const final;
};

/**
 * See (Marson, 2022) pag 101, ELIUL was firstly introduced in (Marson et al., 2021)
 *
 * \Bibliography Marson, F. (2022). Directional lattice Boltzmann boundary conditions.
 * https://doi.org/10.13097/archive-ouverte/unige:160770
 *
 * \Bibliography Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., & Latt, J. (2021).
 * Enhanced single-node lattice Boltzmann boundary condition for fluid flows.
 * Physical Review E, 103(5), 053308. https://doi.org/10.1103/PhysRevE.103.053308
 *
 * @tparam T
 * @tparam D
 */
template <typename T, template <typename U> class D>
class ELIULC : public ELIModels3D<T, D> {
    using ELIModels3D<T, D>::ELIModels3D;
    ELIULC<T, D> *clone() const override;
    inline std::array<T, 5> eliCoefficients(T q, T tauPlus, T tauMinus) const final;
};

/**
 * See (Marson, 2022) pag 102, ELIUL was firstly introduced in (Marson et al., 2021),
 * correction K1 was first introduced in (Ginzburg et al., 2008)
 *
 * Marson, F. (2022). Directional lattice Boltzmann boundary conditions.
 * https://doi.org/10.13097/archive-ouverte/unige:160770
 *
 * Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., & Latt, J. (2021).
 * Enhanced single-node lattice Boltzmann boundary condition for fluid flows.
 * Physical Review E, 103(5), 053308. https://doi.org/10.1103/PhysRevE.103.053308
 *
 * Ginzburg, I., Verhaeghe, F., & d’Humières, D. (2008).
 * Two-Relaxation-Time Lattice Boltzmann Scheme: About Parametrization, Velocity, Pressure and
 * Mixed Boundary Conditions. Commun. Comput. Phys., 3(2), 427–478.
 * @tparam T
 * @tparam D
 */
template <typename T, template <typename U> class D>
class ELIULK1 : public ELIModels3D<T, D> {
    using ELIModels3D<T, D>::ELIModels3D;
    ELIULK1<T, D> *clone() const override;
    inline std::array<T, 5> eliCoefficients(T q, T tauPlus, T tauMinus) const final;
};

/**
 * See (Marson, 2022) pag 102, ELIUL was firstly introduced in (Marson et al., 2021)
 *
 * Marson, F. (2022). Directional lattice Boltzmann boundary conditions.
 * https://doi.org/10.13097/archive-ouverte/unige:160770
 *
 * Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., & Latt, J. (2021).
 * Enhanced single-node lattice Boltzmann boundary condition for fluid flows.
 * Physical Review E, 103(5), 053308. https://doi.org/10.1103/PhysRevE.103.053308
 * @tparam T
 * @tparam D
 */
template <typename T, template <typename U> class D>
class ELIULK3 : public ELIModels3D<T, D> {
    using ELIModels3D<T, D>::ELIModels3D;
    ELIULK3<T, D> *clone() const override;
    inline std::array<T, 5> eliCoefficients(T q, T tauPlus, T tauMinus) const final;
};

/**
 * See (Marson, 2022) pag 102, ELIUL was firstly introduced in (Marson et al., 2021)
 *
 * Marson, F. (2022). Directional lattice Boltzmann boundary conditions.
 * https://doi.org/10.13097/archive-ouverte/unige:160770
 *
 * Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., & Latt, J. (2021).
 * Enhanced single-node lattice Boltzmann boundary condition for fluid flows.
 * Physical Review E, 103(5), 053308. https://doi.org/10.1103/PhysRevE.103.053308
 * @tparam T
 * @tparam D
 */
template <typename T, template <typename U> class D>
class ELIULK4 : public ELIModels3D<T, D> {
    using ELIModels3D<T, D>::ELIModels3D;
    ELIULK4<T, D> *clone() const override;
    inline std::array<T, 5> eliCoefficients(T q, T tauPlus, T tauMinus) const final;
};

// #if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus > 201703L)
// template<typename F>
// concept object = requires(F f){
//                      std::is_object_v<F>(f);
//                  };
// #endif

template <typename T, template <typename U> class D, typename Function>
#if ((defined(_MSVC_LANG) && _MSVC_LANG >= 201703L) || __cplusplus > 201703L)
    requires std::invocable<Function, T, T, T> /*&& object<Function>*/
#endif
class ELIgeneric : public ELIModels3D<T, D> {
public:
    ELIgeneric(BoundaryShape3D<T, Array<T, 3>> *shape_, int flowType_, Function coefficients_);

private:
    ELIgeneric<T, D, Function> *clone() const final;
    inline std::array<T, 5> eliCoefficients(T q, T tauPlus, T tauMinus) const final;
    Function compute_coefficients;
};

}  // namespace plb
#endif  // ELI_OFF_LATTICE_MODEL_3D_H
