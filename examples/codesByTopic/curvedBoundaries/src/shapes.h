/* This code is a showcase for the Palabos library.
 *
 * The Palabos software is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * You can redistribute it and/or modify this code
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 */

#include <math.h>

#include "palabos2D.h"
#include "palabos2D.hh"
#include "simulationParameters.h"

using namespace plb;
using namespace plb::descriptors;
#ifndef CURVED_BOUNDARIES_SHAPES_H
#define CURVED_BOUNDARIES_SHAPES_H

template <typename Real>
TriangleSet<Real> *generateEllipsoid(Array<Real, 3> &center, Real a, Real b,
                                     Real c, Real dx = 1.) {
    Real p = 2.0;  // 1.6075;
    Real approxSurface =
        4. * M_PI *
        pow(1. / 3. * (pow(a * b, p) + pow(a * c, p) + pow(b * c, p)), 1. / p);
    approxSurface /= (dx * dx);
    TriangleSet<Real> sphere =
        constructSphere<Real>(Array<Real, 3>(0, 0, 0), 1., approxSurface);
    sphere.scale(a, b, c);
    sphere.translate(center);
    return sphere.clone();
}
template <typename Real>
TriangleSet<Real> *generateTwoEllipsoid(Array<Real, 3> &center, Real a, Real b,
                                     Real c, Real dx = 1.) {
    Real p = 2.0;
    Real approxSurface =
        4. * M_PI *
        pow(1. / 3. * (pow(a * b, p) + pow(a * c, p) + pow(b * c, p)), 1. / p);
    approxSurface /= (dx * dx);
    TriangleSet<Real> sphere =
        constructSphere<Real>(Array<Real, 3>(0, 0, 0), 1., approxSurface);
    auto& sphere2 = *sphere.clone();
    sphere.scale(a, b, c);
    sphere2.scale(a, b, c);
    sphere.translate(center);
    sphere2.translate(center+Array<Real,3>(5.0*a,.1,.0));
    sphere.append(sphere2);
    return sphere.clone();
}

template <typename Real>
TriangleSet<Real> *readObstacle(Array<Real, 3> &center, Real a, Real b,
                                     Real c, Real scale = 1.) {
    Real p = 2.0;
    TriangleSet<Real> obstacle("curvedBoundaries.stl");
    // Place the obstacle in the correct place in the simulation domain.
    // Here the "geometric center" of the obstacle is computed manually,
    // by computing first its bounding cuboid. In cases that the STL
    // file with the geometry of the obstacle contains its center as
    // the point, say (0, 0, 0), then the following variable
    // "obstacleCenter" must be set to (0, 0, 0) manually.
    Cuboid<Real> bCuboid = obstacle.getBoundingCuboid();
    Array<Real,3> obstacleCenter = (Real) 0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner);
    obstacle.translate(-obstacleCenter);
    obstacle.scale(scale);
    obstacle.rotate(-M_PI/4.,0,0);
    obstacle.scale(a, b, c);
    obstacle.translate(center);
    return obstacle.clone();
}

#endif  // CURVED_BOUNDARIES_SHAPES_H
