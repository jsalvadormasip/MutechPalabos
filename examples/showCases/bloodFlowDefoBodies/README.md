# Instructions for Palabos-npFEM

For more information on the library (numerics & parameters), you can consult the publications below:

* https://arxiv.org/abs/1903.06479
* https://arxiv.org/abs/1911.03062

## Compilation

Palabos-npFEM has been tested in both UNIX and Windows (with the provided CMake File):

1. mkdir build
2. cd build
3. cmake ..
4. make -j && make -j gpu (if ENABLE_CUDA ON)

## Perform Cell Packing

Cell Packing is the software that randomly initializes RBCs & PLTs in the flow field, and after some iterations it resolves all the collisions/interpenetrations. The resolved positions per blood cell are stored in the folder CPs. Open/Explore the cellPacking_params.xml to decide the geometry, hematocrit and many other parameters.

CPU-version:
mpirun -n X ./bloodFlowDefoBodies cellPacking_params.xml

GPU-version:
mpirun -n X ./bloodFlowDefoBodies_gpu cellPacking_params.xml NumberOfNodes NumberOfGPUsPerNode
NumberOfNodes: 1 for a workstation, varies in a cluster
NumberOfGPUsPerNode: number of GPUs per node

In folder tmp you can find the vtks of RBCs & PLTs, along with profiling and other output.

After having the initial positions of the blood cells (stored in the CPs folder), we can perform simulations with a proper flow field.

## Cellular Blood Flow Simuations

These simulations either start from a given CPs folder (generated through Cell Packing or provided) or from a file (initialPlacing.pos) that stores the position (center of mass) and orientation of the blood cells. The initialPlacing.pos stores first the RBCs as ID (0: RBC, 1: PLT), position and orientation (in degrees) and after the PLTs. Cell Packing is the prefered way to generate complex flow fields with many parameters.

CPU-version:
mpirun -n X ./bloodFlowDefoBodies shear_params/poiseuille_params.xml

GPU-version:
mpirun -n X ./bloodFlowDefoBodies_gpu shear_params/poiseuille_params.xml NumberOfNodes NumberOfGPUsPerNode

## Case study: Collision at an obstacle

We provide one case study that uses the initialPlacing.pos. It can be executed as normally with the obstacle_params.xml and it involves 1 RBC colliding with a cubic obstacle.

## Fast Track

Just extract the provided CPs_X.tar.gz and perform a Cellular Blood Flow Simuation as described above. You can find out about the simulation parameters by opening the relevant xml file.
