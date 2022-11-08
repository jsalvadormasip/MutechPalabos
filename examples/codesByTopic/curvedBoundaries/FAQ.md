### Technical FAQ
* **What is `[[maybe_unused]]`** It is a c++17 attribute to suppress unused warning messages from the compiler;
* **I don't have mpi installed** You can run the program in serial, don't use `mpirun -np 2 name_of_exec` but directly `./name_of_exec`
* **What is this return type`auto [lattice_ptr, boundaryon, boundaryoff] = ...`?** It is the c++17 way of unpack elements of a tuple. The tuple is unpacked when the function is returning it.
* **Seg fault after modifying the exercise** Keep in mind that most of the Palabos functions and methods take references as parameters. Thus, if you call them from inside functions, be sure to pass heap allocated variables, or variables allocated in the stack, but inside the main function.
* **How do I tell cmake to build in debug mode?** Use `cmake -DCMAKE_BUILD_TYPE=Debug ..`
* **How do I add a personal stl to the simulation?** We do not cover this topic in this exercise. However, very few modifications are necessary. Please, check the `palabos/examples/showCases/externalFlowAroundObstacle/`
* **Boundary conditions available as June 2021:**
    * Guo scheme (GZS) [1]
    * Bouzidi scheme (BFL) [2]
    * Filippova-Hanel scheme (FH) [3]
    * Mei-Luo-Shyy (MLS) [4]
* **Boundary conditions under validation (unrealeased, available soon):**
    * Central Linear Interpolation scheme (CLI)
    * ELI schemes

### Theoretical FAQ
* **Why don't you show us how to move the boundary?** We do not cover this topic in this exercise for two reasons: for classical boundaries it would be too complicated for 1h exercise (time dependent voxelization, refilling techniques that are not implemented in palabos). For Immersed boundaries there is a good example in Palabos folder.
* **How can have moving boundaries in Palabos?** The easiest way is to use the IBM multi direct forcing approach used in examples/showCases/movingWalls.cpp and in examples/showCases/bloodFlowDefoBodies. If you want to “move” clasical boundary condition the procedure is more complex: 1. Implement with a data processor a refill technique discussed in the theoretical lesson, 2. Voxelize at each time step to track the motion of the surface. Or re-implement everything from scratch using data processors of palabos.
* **What collision model should I use** Depends. If you are working on porous medias and/or low Reynolds flows, it is better to go with a viscosity independent method like CLI and a TRT collision model. At high Reynolds numbers -> RRBGK
* **Which is the most accurate method?** The Multi reflections are the most accurate, but it is not currently implemented in Palabos.
* **What method is used for the force computation** The Momentum Exchange method
* **What single-node methods are implemented in palabos?** For now there is only the FH method. Nevertheless, we are going to release a new version with a variety of single-node methods soon.
