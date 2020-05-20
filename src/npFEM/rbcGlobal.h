#pragma once

namespace plb {

struct LocalMesh
{

    LocalMesh(RawConnectedTriangleMesh<T> const& meshTemplate, pluint bodyID_)
        : mesh(meshTemplate)
        , bodyID(bodyID_)
    { }

    void sendForcesAndCollisionContainers(pluint processorToSend)
    {
        numVertices = (pluint)vertexIDs.size();
        numCollidingNeighbors = (pluint)collisionNeighbors.size();

        // if numVertices == 0: body completely in the envelope.
        // If completely in envelope, no collision neighbors for
        // this local Mesh
        if (numVertices > 0)
        {
#ifdef ENABLE_LOGS
            plb::global::logfile_nonparallel("localProgress.log").flushEntry("MPI_Isend - Sending to : " +
                util::val2str(processorToSend) + ", bodyID=" + util::val2str(bodyID));
#endif // ENABLE_LOGS

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&bodyID, 1, MPI_UNSIGNED_LONG_LONG, (int)processorToSend,
                shapeOpMPItagForcesAndCollisionContainersBodyID, MPI_COMM_WORLD, mpiRequests.back());

#ifdef ENABLE_LOGS
            plb::global::logfile_nonparallel("localProgress.log").flushEntry("Sending numVertices=" + util::val2str(numVertices));
#endif // ENABLE_LOGS
            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&numVertices, 1, MPI_UNSIGNED_LONG_LONG, (int)processorToSend,
                shapeOpMPItagForcesAndCollisionContainersNumVertices, MPI_COMM_WORLD, mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&vertexIDs[0], (int)numVertices, MPI_UNSIGNED_LONG_LONG,
                (int)processorToSend, shapeOpMPItagForcesAndCollisionContainersVertexIDs, MPI_COMM_WORLD,
                mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&shearForces[0][0], (int)(3 * numVertices), MPI_DOUBLE,
                (int)processorToSend, shapeOpMPItagForcesAndCollisionContainersShearForces, MPI_COMM_WORLD,
                mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&normals[0][0], (int)(3 * numVertices), MPI_DOUBLE,
                (int)processorToSend, shapeOpMPItagForcesAndCollisionContainersNormals, MPI_COMM_WORLD,
                mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&pressure[0], (int)numVertices, MPI_DOUBLE,
                (int)processorToSend, shapeOpMPItagForcesAndCollisionContainersPressure, MPI_COMM_WORLD,
                mpiRequests.back());

            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&area[0], (int)numVertices, MPI_DOUBLE,
                (int)processorToSend, shapeOpMPItagForcesAndCollisionContainersArea, MPI_COMM_WORLD,
                mpiRequests.back());


            // Collision Containers
#ifdef ENABLE_LOGS
            plb::global::logfile_nonparallel("localProgress.log").flushEntry("The body has numCollidingNeighbors=" + util::val2str(numCollidingNeighbors));
#endif // ENABLE_LOGS

            // Send this info to avoid bugs and mismatches
            mpiRequests.push_back(new MPI_Request);
            MPI_Isend(&numCollidingNeighbors, 1, MPI_UNSIGNED_LONG_LONG,
                (int)processorToSend, shapeOpMPItagForcesAndCollisionContainersNumCollidingNeighbors, MPI_COMM_WORLD,
                mpiRequests.back());

            if (numCollidingNeighbors > 0)
            {
                mpiRequests.push_back(new MPI_Request);
                MPI_Isend(&collisionNeighbors[0][0],
                    (int)(3 * numCollidingNeighbors), MPI_DOUBLE,
                    (int)processorToSend, shapeOpMPItagForcesAndCollisionContainersCollisionNeighbors, MPI_COMM_WORLD,
                    mpiRequests.back());

                mpiRequests.push_back(new MPI_Request);
                MPI_Isend(&collisionNeighborsNormals[0][0],
                    (int)(3 * numCollidingNeighbors), MPI_DOUBLE,
                    (int)processorToSend, shapeOpMPItagForcesAndCollisionContainerscollisionNeighborsNormals, MPI_COMM_WORLD,
                    mpiRequests.back());
            }
        }
    }

    void completeRequests()
    {
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("bodyID=" + util::val2str(bodyID) +
            ", Start MPI Completion (Forces & Collisions)");
#endif // ENABLE_LOGS

        MPI_Status status;
        for (pluint i = 0; i < (pluint)mpiRequests.size(); ++i) {
            MPI_Wait(mpiRequests[i], &status);
            delete mpiRequests[i];
        }
        mpiRequests.clear();

#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("bodyID=" + util::val2str(bodyID) +
            ", Finish MPI Completion (Forces & Collisions)");
#endif // ENABLE_LOGS
    }

    ///////////////////////////////////////////////////////////////////////////

    // The mesh may refer to a RBC, PLT or any other body
    RawConnectedTriangleMesh<T> mesh; // in lattice units
    pluint bodyID;

    // All fields below refer to the bulk (no envelope)
    // isWrittenTag == true
    pluint numVertices;
    std::vector<pluint> vertexIDs;
    std::vector<Array<T, 3>> shearForces;
    std::vector<Array<T, 3>> normals;
    std::vector<T> pressure;
    std::vector<T> area;

    // Collision Container
    pluint numCollidingNeighbors;
    std::vector<Array<T, 3>> collisionNeighbors, collisionNeighborsNormals;

    // From ShapeOp to Palabos
    std::map<pluint, Array<T, 3>> velocities; // in lattice units

    std::vector<MPI_Request*> mpiRequests; // MPI Completion
};

// Meshes, corresponding to the bodies immersed in the fluid
// in the current processor, used to compute data like area and
// normals. They are reconstructed at every iteration.
// In lattice units.
inline std::map<pluint, LocalMesh*>& LocalMeshes()
{
    // the key of the map is the unique bodyID
    static std::map<pluint, LocalMesh*> instance;
    return instance;
}

}
