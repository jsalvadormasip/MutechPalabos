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

#ifndef IMMERSED_WALLS_2D_HH
#define IMMERSED_WALLS_2D_HH

#include "atomicBlock/dataField2D.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/array.h"
#include "core/globalDefs.h"
#include "immersedWalls2D.h"

namespace plb {

/* ******** ReduceImmersedTorque2D ************************************ */

template <typename T>
ReduceImmersedTorque2D<T>::ReduceImmersedTorque2D(Array<T, 2> const &center_, int reductionFlag_) :
    center(center_),
    sum_torque_ids(Array<plint, 2>(
        this->getStatistics().subscribeSum(), this->getStatistics().subscribeSum())),
    reductionFlag(reductionFlag_)
{ }

template <typename T>
void ReduceImmersedTorque2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);
    Dot2D location = container->getLocation();
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);
    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<Array<T, 2> > const &g = wallData->g;
    std::vector<int> const &flags = wallData->flags;
    PLB_ASSERT(vertices.size() == g.size());
    PLB_ASSERT(vertices.size() == flags.size());

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> vertex = vertices[i];
        if (flags[i] == reductionFlag
            && closedOpenContained(vertex, domain.shift(location.x, location.y))) {
            Array<T, 2> r(vertex - center);
            Array<T, 2> torque(crossProduct(r, g[i]));
            this->getStatistics().gatherSum(sum_torque_ids[0], torque[0]);
            this->getStatistics().gatherSum(sum_torque_ids[1], torque[1]);
        }
    }
}

template <typename T>
ReduceImmersedTorque2D<T> *ReduceImmersedTorque2D<T>::clone() const
{
    return new ReduceImmersedTorque2D<T>(*this);
}

template <typename T>
void ReduceImmersedTorque2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Container Block.
}

template <typename T>
BlockDomain::DomainT ReduceImmersedTorque2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
Array<T, 2> ReduceImmersedTorque2D<T>::getSumTorque() const
{
    return Array<T, 2>(
        this->getStatistics().getSum(sum_torque_ids[0]),
        this->getStatistics().getSum(sum_torque_ids[1]));
}

/* ******** ReduceImmersedForce2D ************************************ */

template <typename T>
ReduceImmersedForce2D<T>::ReduceImmersedForce2D(int reductionFlag_) :
    sum_g_ids(Array<plint, 2>(
        this->getStatistics().subscribeSum(), this->getStatistics().subscribeSum())),
    reductionFlag(reductionFlag_)
{ }

template <typename T>
void ReduceImmersedForce2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);
    Dot2D location = container->getLocation();
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);
    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<Array<T, 2> > const &g = wallData->g;
    std::vector<int> const &flags = wallData->flags;
    PLB_ASSERT(vertices.size() == g.size());
    PLB_ASSERT(vertices.size() == flags.size());

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> vertex = vertices[i];
        if (flags[i] == reductionFlag
            && closedOpenContained(vertex, domain.shift(location.x, location.y))) {
            this->getStatistics().gatherSum(sum_g_ids[0], g[i][0]);
            this->getStatistics().gatherSum(sum_g_ids[1], g[i][1]);
        }
    }
}

template <typename T>
ReduceImmersedForce2D<T> *ReduceImmersedForce2D<T>::clone() const
{
    return new ReduceImmersedForce2D<T>(*this);
}

template <typename T>
void ReduceImmersedForce2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Container Block.
}

template <typename T>
BlockDomain::DomainT ReduceImmersedForce2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
Array<T, 2> ReduceImmersedForce2D<T>::getSumG() const
{
    return Array<T, 2>(
        this->getStatistics().getSum(sum_g_ids[0]), this->getStatistics().getSum(sum_g_ids[1]));
}

/* ******** ReduceImmersedArea2D ************************************ */

template <typename T>
ReduceImmersedArea2D<T>::ReduceImmersedArea2D(int reductionFlag_) :
    sum_area_id(this->getStatistics().subscribeSum()), reductionFlag(reductionFlag_)
{ }

template <typename T>
void ReduceImmersedArea2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);
    Dot2D location = container->getLocation();
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);
    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    std::vector<int> const &flags = wallData->flags;
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(vertices.size() == flags.size());

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> vertex = vertices[i];
        if (flags[i] == reductionFlag
            && closedOpenContained(vertex, domain.shift(location.x, location.y))) {
            this->getStatistics().gatherSum(sum_area_id, areas[i]);
        }
    }
}

template <typename T>
ReduceImmersedArea2D<T> *ReduceImmersedArea2D<T>::clone() const
{
    return new ReduceImmersedArea2D<T>(*this);
}

template <typename T>
void ReduceImmersedArea2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Container Block.
}

template <typename T>
BlockDomain::DomainT ReduceImmersedArea2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

template <typename T>
T ReduceImmersedArea2D<T>::getSumArea() const
{
    return this->getStatistics().getSum(sum_area_id);
}

/* ******** InamuroIteration2D ************************************ */

template <typename T, class VelFunction>
InamuroIteration2D<T, VelFunction>::InamuroIteration2D(
    VelFunction velFunction_, T tau_, bool incompressibleModel_) :
    velFunction(velFunction_), tau(tau_), incompressibleModel(incompressibleModel_)
{ }

template <typename T, class VelFunction>
void InamuroIteration2D<T, VelFunction>::processGenericBlocks(
    Box2D, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[0]);
    TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[1]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[2]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);
    PLB_ASSERT(container);
    Dot2D location = rhoBar->getLocation();
    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<Array<T, 2> > deltaG(vertices.size());
    std::vector<Array<T, 2> > &g = wallData->g;
    PLB_ASSERT(vertices.size() == g.size());

    // In this iteration, the force is computed for every vertex.
    if (incompressibleModel) {
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            // Use the weighting function to compute the average momentum
            // and the average density on the surface vertex.
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(vertex);
            deltaG[i] = areas[i] * (wallVelocity - averageJ);
            g[i] += deltaG[i];
        }
    } else {  // Compressible model.
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            T averageRhoBar = T();
            // Use the weighting function to compute the average momentum
            // and the average density on the surface vertex.
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                    averageRhoBar += W * nextRhoBar;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(vertex);
            deltaG[i] = areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
            // g[i] += deltaG[i];
            g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
        }
    }

    // In this iteration, the force is applied from every vertex to the grid nodes.
    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
            for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                nextJ += tau * W * deltaG[i];
                j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
            }
        }
    }
}

template <typename T, class VelFunction>
InamuroIteration2D<T, VelFunction> *InamuroIteration2D<T, VelFunction>::clone() const
{
    return new InamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void InamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT InamuroIteration2D<T, VelFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** IndexedInamuroIteration2D ************************************ */

template <typename T, class VelFunction>
IndexedInamuroIteration2D<T, VelFunction>::IndexedInamuroIteration2D(
    VelFunction velFunction_, T tau_, bool incompressibleModel_) :
    velFunction(velFunction_), tau(tau_), incompressibleModel(incompressibleModel_)
{ }

template <typename T, class VelFunction>
void IndexedInamuroIteration2D<T, VelFunction>::processGenericBlocks(
    Box2D, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[0]);
    TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[1]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[2]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);
    PLB_ASSERT(container);
    Dot2D location = rhoBar->getLocation();
    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<Array<T, 2> > deltaG(vertices.size());
    std::vector<Array<T, 2> > &g = wallData->g;
    PLB_ASSERT(vertices.size() == g.size());
    std::vector<pluint> const &globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT(vertices.size() == globalVertexIds.size());

    if (incompressibleModel) {
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(globalVertexIds[i]);
            deltaG[i] = areas[i] * (wallVelocity - averageJ);
            g[i] += deltaG[i];
        }
    } else {  // Compressible model.
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            T averageRhoBar = T();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
                for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    if (nextRhoBar == -1)
                        continue;
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                    averageRhoBar += W * nextRhoBar;
                }
            }
            // averageJ += (T)0.5*g[i];
            Array<T, 2> wallVelocity = velFunction(globalVertexIds[i]);
            deltaG[i] = areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
            // g[i] += deltaG[i];
            g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
        }
    }

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        for (plint dx = xLim[0]; dx <= xLim[1]; ++dx) {
            for (plint dy = yLim[0]; dy <= yLim[1]; ++dy) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                //	WARN: May be wrong! There must be complete grid point in the range of  delta
                // function!
                //                T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                //                if (nextRhoBar == -1)
                //                    continue;

                Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                nextJ += tau * W * deltaG[i];
                j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
            }
        }
    }
}

template <typename T, class VelFunction>
IndexedInamuroIteration2D<T, VelFunction> *IndexedInamuroIteration2D<T, VelFunction>::clone() const
{
    return new IndexedInamuroIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void IndexedInamuroIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT IndexedInamuroIteration2D<T, VelFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** IndexedImmersedBoundaryIteration2D ************************************ */

template <typename T, class VelFunction>
IndexedImmersedBoundaryIteration2D<T, VelFunction>::IndexedImmersedBoundaryIteration2D(
    VelFunction velFunction_) :
    velFunction(velFunction_)
{ }

template <typename T, class VelFunction>
void IndexedImmersedBoundaryIteration2D<T, VelFunction>::processGenericBlocks(
    [[maybe_unused]] Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    TensorField2D<T, 2> *u = dynamic_cast<TensorField2D<T, 2> *>(blocks[0]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[1]);
    PLB_ASSERT(u);
    PLB_ASSERT(container);
    Dot2D location = u->getLocation();
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<Array<T, 2> > deltaG(vertices.size());
    std::vector<Array<T, 2> > &g = wallData->g;
    PLB_ASSERT(vertices.size() == g.size());
    std::vector<pluint> const &globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT(vertices.size() == globalVertexIds.size());

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        Array<T, 2> averageU;
        averageU.resetToZero();
        // x   x . x   x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                Array<T, 2> nextU = u->get(pos[0], pos[1]);
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                averageU += W * nextU;
            }
        }
        // averageU += (T)0.5*g[i];
        Array<T, 2> wallVelocity = velFunction(globalVertexIds[i]);
        deltaG[i] = areas[i] * (wallVelocity - averageU);
        g[i] += (T)2 * deltaG[i];
    }

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                Array<T, 2> nextU = u->get(pos[0], pos[1]);
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                nextU += W * deltaG[i];
                u->get(pos[0], pos[1]) = nextU;
            }
        }
    }
}

template <typename T, class VelFunction>
IndexedImmersedBoundaryIteration2D<T, VelFunction>
    *IndexedImmersedBoundaryIteration2D<T, VelFunction>::clone() const
{
    return new IndexedImmersedBoundaryIteration2D<T, VelFunction>(*this);
}

template <typename T, class VelFunction>
void IndexedImmersedBoundaryIteration2D<T, VelFunction>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Velocity
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, class VelFunction>
BlockDomain::DomainT IndexedImmersedBoundaryIteration2D<T, VelFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ConstVelInamuroIteration2D ************************************ */

template <typename T>
ConstVelInamuroIteration2D<T>::ConstVelInamuroIteration2D(
    Array<T, 2> const &wallVelocity_, T tau_, bool incompressibleModel_) :
    wallVelocity(wallVelocity_), tau(tau_), incompressibleModel(incompressibleModel_)
{ }

template <typename T>
void ConstVelInamuroIteration2D<T>::processGenericBlocks(
    [[maybe_unused]] Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);
    ScalarField2D<T> *rhoBar = dynamic_cast<ScalarField2D<T> *>(blocks[0]);
    TensorField2D<T, 2> *j = dynamic_cast<TensorField2D<T, 2> *>(blocks[1]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[2]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(j);
    PLB_ASSERT(container);
    Dot2D location = rhoBar->getLocation();
    Dot2D ofsJ = computeRelativeDisplacement(*rhoBar, *j);
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);
    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<T> const &areas = wallData->areas;
    PLB_ASSERT(vertices.size() == areas.size());
    std::vector<Array<T, 2> > deltaG(vertices.size());
    std::vector<Array<T, 2> > &g = wallData->g;
    PLB_ASSERT(vertices.size() == g.size());

    if (incompressibleModel) {
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                }
            }
            // averageJ += (T)0.5*g[i];
            deltaG[i] = areas[i] * (wallVelocity - averageJ);
            g[i] += deltaG[i];
        }
    } else {  // Compressible model.
        for (pluint i = 0; i < vertices.size(); ++i) {
            Array<T, 2> const &vertex = vertices[i];
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
            Array<T, 2> averageJ;
            averageJ.resetToZero();
            T averageRhoBar = T();
            // x   x . x   x
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1]);
                    Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                    Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                    T W = inamuroDeltaFunction2D<T>().W(r);
                    averageJ += W * nextJ;
                    averageRhoBar += W * nextRhoBar;
                }
            }
            // averageJ += (T)0.5*g[i];
            deltaG[i] = areas[i] * ((averageRhoBar + (T)1.) * wallVelocity - averageJ);
            // g[i] += deltaG[i];
            g[i] += deltaG[i] / ((T)1.0 + averageRhoBar);
        }
    }

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                Array<T, 2> nextJ = j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y);
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                nextJ += tau * W * deltaG[i];
                j->get(pos[0] + ofsJ.x, pos[1] + ofsJ.y) = nextJ;
            }
        }
    }
}

template <typename T>
ConstVelInamuroIteration2D<T> *ConstVelInamuroIteration2D<T>::clone() const
{
    return new ConstVelInamuroIteration2D<T>(*this);
}

template <typename T>
void ConstVelInamuroIteration2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT ConstVelInamuroIteration2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ComputeImmersedBoundaryForce2D ************************************ */

template <typename T>
void ComputeImmersedBoundaryForce2D<T>::processGenericBlocks(
    [[maybe_unused]] Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    TensorField2D<T, 2> *force = dynamic_cast<TensorField2D<T, 2> *>(blocks[0]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[1]);
    PLB_ASSERT(force);
    PLB_ASSERT(container);
    Dot2D location = force->getLocation();
    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);
    std::vector<Array<T, 2> > const &vertices = wallData->vertices;
    std::vector<Array<T, 2> > &g = wallData->g;
    PLB_ASSERT(vertices.size() == g.size());

    for (plint iX = 0; iX < force->getNx(); iX++) {
        for (plint iY = 0; iY < force->getNy(); iY++) {
            force->get(iX, iY).resetToZero();
        }
    }

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
        const Array<T, 2> fracPos(util::frac(vertex[0]), util::frac(vertex[1]));
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                Array<T, 2> r((T)dx - fracPos[0], (T)dy - fracPos[1]);
                T W = inamuroDeltaFunction2D<T>().W(r);
                force->get(pos[0], pos[1]) += W * g[i];
            }
        }
    }
}

template <typename T>
ComputeImmersedBoundaryForce2D<T> *ComputeImmersedBoundaryForce2D<T>::clone() const
{
    return new ComputeImmersedBoundaryForce2D<T>(*this);
}

template <typename T>
void ComputeImmersedBoundaryForce2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Force
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT ComputeImmersedBoundaryForce2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallData2D ************************************ */

template <typename T>
InstantiateImmersedWallData2D<T>::InstantiateImmersedWallData2D(
    std::vector<Array<T, 2> > const &vertices_, std::vector<T> const &areas_,
    std::vector<Array<T, 2> > const &normals_) :
    vertices(vertices_), areas(areas_), normals(normals_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(normals.size() == 0 || normals.size() == areas.size());
}

template <typename T>
void InstantiateImmersedWallData2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);
    bool useNormals = normals.size() > 0;
    Dot2D location = container->getLocation();
    ImmersedWallData2D<T> *wallData = new ImmersedWallData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(2).shift(location.x, location.y));

    for (pluint i = 0; i < vertices.size(); ++i) {
        // Vertices which are close to the boundaries of the extendedEnvelope
        // are irrelevant, because they will act upon the bulk of the computational
        // domain through an Inamuro kernel, which at this distance is close to zero.
        // It is therefore OK, numerically speaking to exclude an epsilon-margin close
        // to these boundaries. Plus, it is required for technical reasons, because if
        // later on we pass across the boundaries of the extendedEnvelope because
        // of roundoff errors, the code will crash.
        static const T epsilon = 1.e-4;
        if (contained(vertices[i], extendedEnvelope, epsilon)) {
            wallData->vertices.push_back(vertices[i]);
            wallData->areas.push_back(areas[i]);
            if (useNormals) {
                wallData->normals.push_back(normals[i]);
            }
            wallData->g.push_back(Array<T, 2>((T)0., (T)0.));
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->flags = std::vector<int>(wallData->vertices.size(), 0);
    container->setData(wallData);
}

template <typename T>
InstantiateImmersedWallData2D<T> *InstantiateImmersedWallData2D<T>::clone() const
{
    return new InstantiateImmersedWallData2D<T>(*this);
}

template <typename T>
void InstantiateImmersedWallData2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedWallData2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallDataWithTagging2D ************************************ */

template <typename T>
InstantiateImmersedWallDataWithTagging2D<T>::InstantiateImmersedWallDataWithTagging2D(
    std::vector<Array<T, 2> > const &vertices_, std::vector<T> const &areas_, int fluidFlag_) :
    vertices(vertices_), areas(areas_), fluidFlag(fluidFlag_)
{
    PLB_ASSERT(vertices.size() == areas.size());
}

template <typename T>
void InstantiateImmersedWallDataWithTagging2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);
    Dot2D location = container->getLocation();

    ScalarField2D<int> *flagMatrix = dynamic_cast<ScalarField2D<int> *>(blocks[1]);
    PLB_ASSERT(flagMatrix);
    Dot2D ofsFlag = computeRelativeDisplacement(*container, *flagMatrix);
    Array<plint, 2> flagDispl(ofsFlag.x, ofsFlag.y);

    ImmersedWallData2D<T> *wallData = new ImmersedWallData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(2).shift(location.x, location.y));

    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> vertex = vertices[i];
        // Vertices which are close to the boundaries of the extendedEnvelope
        // are irrelevant, because they will act upon the bulk of the computational
        // domain through an Inamuro kernel, which at this distance is close to zero.
        // It is therefore OK, numerically speaking to exclude an epsilon-margin close
        // to these boundaries. Plus, it is required for technical reasons, because if
        // later on we pass across the boundaries of the extendedEnvelope because
        // of roundoff errors, the code will crash.
        static const T epsilon = 1.e-4;
        if (contained(vertex, extendedEnvelope, epsilon)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            wallData->g.push_back(Array<T, 2>((T)0., (T)0.));
            wallData->globalVertexIds.push_back(i);
            Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
            const Array<plint, 2> xLim(
                (vertex[0] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            const Array<plint, 2> yLim(
                (vertex[1] < (T)0 ? Array<plint, 2>(-2, 1) : Array<plint, 2>(-1, 2)));
            bool hasFluidNeighbor = false;
            for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
                for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                    Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy) + flagDispl);
                    if (flagMatrix->get(pos[0], pos[1]) == fluidFlag) {
                        hasFluidNeighbor = true;
                    }
                }
            }
            if (hasFluidNeighbor) {
                wallData->flags.push_back(0);
            } else {
                wallData->flags.push_back(1);
            }
        }
    }
    container->setData(wallData);
}

template <typename T>
InstantiateImmersedWallDataWithTagging2D<T> *InstantiateImmersedWallDataWithTagging2D<T>::clone()
    const
{
    return new InstantiateImmersedWallDataWithTagging2D<T>(*this);
}

template <typename T>
void InstantiateImmersedWallDataWithTagging2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
    modified[1] = modif::nothing;          // Flag matrix.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedWallDataWithTagging2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}
/* ******** InstantiateImmersedWallDataWithIndexedTagging2D ************************************ */

template <typename T>
InstantiateImmersedWallDataWithIndexedTagging2D<T>::InstantiateImmersedWallDataWithIndexedTagging2D(
    std::vector<Array<T, 2> > const &vertices_, std::vector<T> const &areas_,
    std::vector<int> const &flags_) :
    vertices(vertices_), areas(areas_), flags(flags_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(vertices.size() == flags.size());
}

template <typename T>
void InstantiateImmersedWallDataWithIndexedTagging2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);
    Dot2D location = container->getLocation();

    ImmersedWallData2D<T> *wallData = new ImmersedWallData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(2).shift(location.x, location.y));

    for (pluint i = 0; i < vertices.size(); ++i) {
        // Vertices which are close to the boundaries of the extendedEnvelope
        // are irrelevant, because they will act upon the bulk of the computational
        // domain through an Inamuro kernel, which at this distance is close to zero.
        // It is therefore OK, numerically speaking to exclude an epsilon-margin close
        // to these boundaries. Plus, it is required for technical reasons, because if
        // later on we pass across the boundaries of the extendedEnvelope because
        // of roundoff errors, the code will crash.
        static const T epsilon = 1.e-4;
        if (contained(vertices[i], extendedEnvelope, epsilon)) {
            wallData->vertices.push_back(vertices[i]);
            wallData->areas.push_back(areas[i]);
            wallData->g.push_back(Array<T, 2>((T)0., (T)0.));
            wallData->flags.push_back(flags[i]);
            wallData->globalVertexIds.push_back(i);
        }
    }
    container->setData(wallData);
}

template <typename T>
InstantiateImmersedWallDataWithIndexedTagging2D<T>
    *InstantiateImmersedWallDataWithIndexedTagging2D<T>::clone() const
{
    return new InstantiateImmersedWallDataWithIndexedTagging2D<T>(*this);
}

template <typename T>
void InstantiateImmersedWallDataWithIndexedTagging2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT InstantiateImmersedWallDataWithIndexedTagging2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** InstantiateSurfaceBlockData2D ************************************ */

template <typename T>
InstantiateSurfaceBlockData2D<T>::InstantiateSurfaceBlockData2D(
    plint envelopeWidth_, std::vector<Array<T, 2> > const &vertices_) :
    envelopeWidth(envelopeWidth_), vertices(vertices_)
{
    PLB_ASSERT(envelopeWidth >= 1);
}

template <typename T>
void InstantiateSurfaceBlockData2D<T>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);

    Dot2D location = container->getLocation();

    SurfaceBlockData2D<T> *surfaceData = new SurfaceBlockData2D<T>;
    Box2D extendedEnvelope(domain.enlarge(envelopeWidth).shift(location.x, location.y));

    for (pluint i = 0; i < vertices.size(); ++i) {
        if (containedInclusive(vertices[i], extendedEnvelope)) {
            surfaceData->vertices.push_back(vertices[i]);
        }
    }
    container->setData(surfaceData);
}

template <typename T>
InstantiateSurfaceBlockData2D<T> *InstantiateSurfaceBlockData2D<T>::clone() const
{
    return new InstantiateSurfaceBlockData2D<T>(*this);
}

template <typename T>
void InstantiateSurfaceBlockData2D<T>::getTypeOfModification(
    std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT InstantiateSurfaceBlockData2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** SurfaceOnLattice2D *************************************** */

template <typename T, typename U>
SurfaceOnLattice2D<T, U>::SurfaceOnLattice2D(U value_, plint envelopeWidth_) :
    value(value_), envelopeWidth(envelopeWidth_)
{
    PLB_ASSERT(envelopeWidth >= 1);
}

template <typename T, typename U>
void SurfaceOnLattice2D<T, U>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    ScalarField2D<U> *surfaceOnLattice = dynamic_cast<ScalarField2D<U> *>(blocks[0]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[1]);
    PLB_ASSERT(surfaceOnLattice);
    PLB_ASSERT(container);

    SurfaceBlockData2D<T> *surfaceData =
        dynamic_cast<SurfaceBlockData2D<T> *>(container->getData());
    PLB_ASSERT(surfaceData);

    Dot2D location = surfaceOnLattice->getLocation();

    std::vector<Array<T, 2> > const &vertices = surfaceData->vertices;

    plint w = envelopeWidth;
    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-(w - 1) - 1, w - 1)
                              : Array<plint, 2>(-(w - 1), w)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-(w - 1) - 1, w - 1)
                              : Array<plint, 2>(-(w - 1), w)));
        // x  x  x . x  x  x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                if (contained(pos[0], pos[1], domain)) {
                    surfaceOnLattice->get(pos[0], pos[1]) = value;
                }
            }
        }
    }
}

template <typename T, typename U>
SurfaceOnLattice2D<T, U> *SurfaceOnLattice2D<T, U>::clone() const
{
    return new SurfaceOnLattice2D<T, U>(*this);
}

template <typename T, typename U>
void SurfaceOnLattice2D<T, U>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Surface values (ScalarField2D<U>).
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, typename U>
BlockDomain::DomainT SurfaceOnLattice2D<T, U>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** SurfaceOnLattice2D_N *************************************** */

template <typename T, typename U>
SurfaceOnLattice2D_N<T, U>::SurfaceOnLattice2D_N(U value_, plint envelopeWidth_) :
    value(value_), envelopeWidth(envelopeWidth_)
{
    PLB_ASSERT(envelopeWidth >= 1);
}

template <typename T, typename U>
void SurfaceOnLattice2D_N<T, U>::processGenericBlocks(
    Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);
    NTensorField2D<U> *surfaceOnLattice = dynamic_cast<NTensorField2D<U> *>(blocks[0]);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[1]);
    PLB_ASSERT(surfaceOnLattice);
    PLB_ASSERT(container);

    SurfaceBlockData2D<T> *surfaceData =
        dynamic_cast<SurfaceBlockData2D<T> *>(container->getData());
    PLB_ASSERT(surfaceData);

    Dot2D location = surfaceOnLattice->getLocation();

    std::vector<Array<T, 2> > const &vertices = surfaceData->vertices;

    plint w = envelopeWidth;
    for (pluint i = 0; i < vertices.size(); ++i) {
        Array<T, 2> const &vertex = vertices[i];
        Array<plint, 2> intPos((plint)vertex[0] - location.x, (plint)vertex[1] - location.y);
        const Array<plint, 2> xLim(
            (vertex[0] < (T)0 ? Array<plint, 2>(-(w - 1) - 1, w - 1)
                              : Array<plint, 2>(-(w - 1), w)));
        const Array<plint, 2> yLim(
            (vertex[1] < (T)0 ? Array<plint, 2>(-(w - 1) - 1, w - 1)
                              : Array<plint, 2>(-(w - 1), w)));
        // x  x  x . x  x  x
        for (plint dx = xLim[0]; dx <= xLim[1]; dx++) {
            for (plint dy = yLim[0]; dy <= yLim[1]; dy++) {
                Array<plint, 2> pos(intPos + Array<plint, 2>(dx, dy));
                if (contained(pos[0], pos[1], domain)) {
                    *surfaceOnLattice->get(pos[0], pos[1]) = value;
                }
            }
        }
    }
}

template <typename T, typename U>
SurfaceOnLattice2D_N<T, U> *SurfaceOnLattice2D_N<T, U>::clone() const
{
    return new SurfaceOnLattice2D_N<T, U>(*this);
}

template <typename T, typename U>
void SurfaceOnLattice2D_N<T, U>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::staticVariables;  // Surface values (NTensorField3D<U>).
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template <typename T, typename U>
BlockDomain::DomainT SurfaceOnLattice2D_N<T, U>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** ResetForceStatistics2D ************************************ */

template <typename T>
void ResetForceStatistics2D<T>::processGenericBlocks(
    [[maybe_unused]] Box2D domain, std::vector<AtomicBlock2D *> blocks)
{
    PLB_PRECONDITION(blocks.size() == 1);
    AtomicContainerBlock2D *container = dynamic_cast<AtomicContainerBlock2D *>(blocks[0]);
    PLB_ASSERT(container);

    ImmersedWallData2D<T> *wallData = dynamic_cast<ImmersedWallData2D<T> *>(container->getData());
    PLB_ASSERT(wallData);

    std::vector<Array<T, 2> > &g = wallData->g;

    for (pluint i = 0; i < g.size(); i++) {
        g[i].resetToZero();
    }
}

template <typename T>
ResetForceStatistics2D<T> *ResetForceStatistics2D<T>::clone() const
{
    return new ResetForceStatistics2D<T>(*this);
}

template <typename T>
void ResetForceStatistics2D<T>::getTypeOfModification(std::vector<modif::ModifT> &modified) const
{
    modified[0] = modif::nothing;  // Container Block with triangle data.
}

template <typename T>
BlockDomain::DomainT ResetForceStatistics2D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}
}  // namespace plb

#endif  // IMMERSED_WALLS_2D_HH
