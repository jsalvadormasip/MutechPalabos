#pragma once

namespace plb {
namespace npfem {

// Do not use MPI_TAG = 0 because it is reserved for Palabos
const int MPItagIniPalabosParticles = 2709;

const int shapeOpMPItagForcesAndCollisionContainersBodyID = 2710;
const int shapeOpMPItagForcesAndCollisionContainersNumVertices = 2711;
const int shapeOpMPItagForcesAndCollisionContainersVertexIDs = 2712;
const int shapeOpMPItagForcesAndCollisionContainersShearForces = 2713;
const int shapeOpMPItagForcesAndCollisionContainersNormals = 2714;
const int shapeOpMPItagForcesAndCollisionContainersPressure = 2715;
const int shapeOpMPItagForcesAndCollisionContainersArea = 2716;
const int shapeOpMPItagForcesAndCollisionContainersNumCollidingNeighbors = 2717;
const int shapeOpMPItagForcesAndCollisionContainersCollisionNeighbors = 2718;
const int shapeOpMPItagForcesAndCollisionContainerscollisionNeighborsNormals = 2719;

const int shapeOpMPItagVelocities = 2730;

}
}