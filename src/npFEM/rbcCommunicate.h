#pragma once

namespace plb
{

// This is the "send" part of the MPI communication for forces (Palabos->ShapeOp) & Collision Containers
template <typename T>
void sendBodyForcesAndCollisionData(std::map<pluint, pluint>& bodyToProc, std::vector<ShapeOpBody>& shapeOpBodies)
{
    for (typename std::map<pluint, LocalMesh*>::iterator it = LocalMeshes().begin(); it != LocalMeshes().end(); ++it)
    {
        it->second->numVertices = (pluint)it->second->vertexIDs.size();

#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log")
            .flushEntry("Local Mesh with ID=" + util::val2str(it->first) +
                " | numVertices=" + util::val2str(it->second->numVertices) + ", at Processor=" + util::val2str(global::mpi().getRank()) +
                ", to Processor=" + util::val2str(bodyToProc[it->second->bodyID]));
#endif // ENABLE_LOGS

        // numVertices == 0 means that the LocalMesh is
        // completely in the envelope
        // No need for anything
        if (it->second->numVertices == 0)
            continue;

        // The bodies on the ShapeOp side are distributed once at the beginning
        // and stay there for the whole execution of the program.
        // Thus, we know at any time where to send the containers.
        pluint processorToSend = bodyToProc[it->second->bodyID];

        if (processorToSend == (pluint)global::mpi().getRank())
        {
            pluint solverID = getSolverID(it->second->bodyID, shapeOpBodies);
            
            shapeOpBodies[solverID].receiveContainers(it->second->vertexIDs,
                it->second->shearForces, it->second->normals,
                it->second->pressure, it->second->area, processorToSend,
                it->second->collisionNeighbors, it->second->collisionNeighborsNormals);
        }
        else
        {   
            it->second->sendForcesAndCollisionContainers(processorToSend);
        }
    }
}

// This is the "receive" part of the MPI communication for forces & Collisions
// Containers (Palabos->ShapeOp).
template <typename T>
void receiveBodyForcesAndCollisionData(std::vector<ShapeOpBody>& shapeOpBodies)
{
    pluint numLocalBodies = shapeOpBodies.size();

    bool isComplete;
    if (numLocalBodies == 0)
        return;
    else
        isComplete = communicationComplete(shapeOpBodies);

#ifdef ENABLE_LOGS
    plb::global::logfile_nonparallel("localProgress.log").flushEntry("\n");
    plb::global::logfile_nonparallel("localProgress.log").flushEntry("RECV Forces & Collision Containers");
    if (isComplete)
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Communication Complete - NO MPI");
    else
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Communication NOT Complete - MPI");
#endif // ENABLE_LOGS

    while (!isComplete)
    {
        pluint bodyID, numVertices;
        MPI_Status status, dummyStatus;

#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("MPI_Recv - MPI_ANY_SOURCE");
#endif // ENABLE_LOGS

        MPI_Recv(&bodyID, 1, MPI_UNSIGNED_LONG_LONG, MPI_ANY_SOURCE,
            shapeOpMPItagForcesAndCollisionContainersBodyID, MPI_COMM_WORLD, &status);
        int sending_processor = status.MPI_SOURCE;

#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("MPI_Recv - From: " +
            util::val2str(sending_processor) + ", bodyID=" + util::val2str(bodyID));
#endif // ENABLE_LOGS

        pluint solverID = getSolverID(bodyID, shapeOpBodies);

        MPI_Recv(&numVertices, 1, MPI_UNSIGNED_LONG_LONG, sending_processor,
            shapeOpMPItagForcesAndCollisionContainersNumVertices, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("numVertices=" + util::val2str(numVertices));
#endif // ENABLE_LOGS

        std::vector<pluint> vertexIDs(numVertices);
        std::vector<Array<T, 3>> shearForces(numVertices);
        std::vector<Array<T, 3>> normals(numVertices);
        std::vector<T> pressure(numVertices);
        std::vector<T> area(numVertices);

        MPI_Recv(&vertexIDs[0], (int)numVertices, MPI_UNSIGNED_LONG_LONG,
            sending_processor, shapeOpMPItagForcesAndCollisionContainersVertexIDs, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received vertexID");
#endif // ENABLE_LOGS

        MPI_Recv(&shearForces[0][0], (int)(3 * numVertices), MPI_DOUBLE,
            sending_processor, shapeOpMPItagForcesAndCollisionContainersShearForces, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received shearForces");
#endif // ENABLE_LOGS

        MPI_Recv(&normals[0][0], (int)(3 * numVertices), MPI_DOUBLE,
            sending_processor, shapeOpMPItagForcesAndCollisionContainersNormals, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received normals");
#endif // ENABLE_LOGS

        MPI_Recv(&pressure[0], (int)(numVertices), MPI_DOUBLE,
            sending_processor, shapeOpMPItagForcesAndCollisionContainersPressure, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received pressure");
#endif // ENABLE_LOGS

        MPI_Recv(&area[0], (int)(numVertices), MPI_DOUBLE, sending_processor,
            shapeOpMPItagForcesAndCollisionContainersArea, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received area");
#endif // ENABLE_LOGS

        
        // Collision Containers (if any)

        
        pluint numCollidingNeighbors;
        MPI_Recv(&numCollidingNeighbors, 1, MPI_UNSIGNED_LONG_LONG,
            sending_processor, shapeOpMPItagForcesAndCollisionContainersNumCollidingNeighbors, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
        plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received numCollidingNeighbors=" +
        util::val2str(numCollidingNeighbors));
#endif // ENABLE_LOGS

        std::vector<Array<T, 3>> collisionNeighbors(numCollidingNeighbors);
        std::vector<Array<T, 3>> collisionNeighborsNormals(numCollidingNeighbors);
        if (numCollidingNeighbors > 0)
        {
            MPI_Recv(&collisionNeighbors[0][0],
                (int)(3 * numCollidingNeighbors), MPI_DOUBLE, sending_processor,
                shapeOpMPItagForcesAndCollisionContainersCollisionNeighbors, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
            plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received neighbors");
#endif // ENABLE_LOGS

            MPI_Recv(&collisionNeighborsNormals[0][0],
                (int)(3 * numCollidingNeighbors), MPI_DOUBLE, sending_processor,
                shapeOpMPItagForcesAndCollisionContainerscollisionNeighborsNormals, MPI_COMM_WORLD, &dummyStatus);
#ifdef ENABLE_LOGS
            plb::global::logfile_nonparallel("localProgress.log").flushEntry("Received neighborsNormals");
#endif // ENABLE_LOGS
        }


        shapeOpBodies[solverID].receiveContainers(vertexIDs, shearForces, normals, pressure, area, sending_processor,
            collisionNeighbors, collisionNeighborsNormals);

        
        // Stop receiving ?
        isComplete = communicationComplete(shapeOpBodies);
#ifdef ENABLE_LOGS
        if (isComplete)
            plb::global::logfile_nonparallel("localProgress.log").flushEntry("Communication Complete");
        else
            plb::global::logfile_nonparallel("localProgress.log").flushEntry("Communication NOT Complete - Re-enter MPI Loop");
#endif // ENABLE_LOGS
    }

#ifdef ENABLE_LOGS
    plb::global::logfile_nonparallel("localProgress.log").flushEntry("\n");
#endif // ENABLE_LOGS
}

// This is the "send" part of the MPI communication for velocities
// (ShapeOp->Palabos).
template <typename T>
void sendBodiesVelocities(std::vector<ShapeOpBody>& shapeOpBodies)
{    
    for (pluint i = 0; i < shapeOpBodies.size(); ++i)
        shapeOpBodies[i].sendVelocities();
}

// This is the "receive" part of the MPI communication for velocities (ShapeOp->Palabos).
template <typename T>
void receiveBodiesVelocities(std::map<pluint, pluint>& bodyToProc, std::vector<ShapeOpBody>& shapeOpBodies, T dx, T dt, T rho)
{
    // From physical to lattice units (From ShapeOp to Palabos)
    T Cu = dt / dx;

    typename std::map<pluint, LocalMesh*>::iterator it = LocalMeshes().begin();
    for (; it != LocalMeshes().end(); ++it)
    {
        if (it->second->numVertices == 0)
            continue;

        pluint bodyID = it->first;              
        int sendingProcessor = bodyToProc[bodyID];
        
        // prepare containers (merely)
        std::vector<Array<T, 3>> velocities;
        it->second->velocities.clear();

        if (sendingProcessor == global::mpi().getRank())
        {
            pluint solverID = getSolverID(bodyID, shapeOpBodies);
            for (pluint j = 0; j < shapeOpBodies[solverID].numTasksToSendVelocities(); ++j) {
                if (shapeOpBodies[solverID].getTask(j) == (pluint)global::mpi().getRank()) {
                    shapeOpBodies[solverID].fillVelocityContainers(j, velocities);
                    break;
                }
            }
        } 
        else // Treating MPI Case
        {
            MPI_Status dummyStatus;

            // prepare containers
            velocities.clear();
            velocities.resize(it->second->numVertices);

            MPI_Recv(&velocities[0][0], (int)(3 * it->second->numVertices), MPI_DOUBLE,
                sendingProcessor, shapeOpMPItagVelocities, MPI_COMM_WORLD, &dummyStatus);
        }

        
        for (pluint i = 0; i < (pluint)it->second->vertexIDs.size(); ++i)
        {
            // Each element of vertexIDs corresponds to each elem of velocities
            // one-to-one correspondence
            pluint vertexID = it->second->vertexIDs[i];
            it->second->velocities[vertexID] = Array<T, 3>(velocities[i][0], velocities[i][1], velocities[i][2]) * Cu;
        }
    }
}

}
