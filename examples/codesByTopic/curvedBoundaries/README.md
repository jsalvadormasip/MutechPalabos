---
title: "Link-wise (directional) boundary conditions for LBM"
author: Francesco Marson, Department of Computer Science, University of Geneva
date: Palabos Summer School 2021
geometry: margin=2cm
output: pdf_document
---

The main goal of this example is to show how you can try and test different curved boundary conditions.
In particular, it shows how to use the the genericELI class that allows to define any LI+ methods.

LI+ methods includes two types of boundaries:
* LI (e.g. Bouzidi et al. 2001, or Yu et al. 2003, or CLI Ginzburg et al. 2008)
* ELI: a class of local methods introduced in (Marson et al. 2021) and extended in (Marson 2022 and Ginzburg et al. 2022)

As beneficial collateral effects, you may familiarize yourself with the following topics

* Compile an example with CMake;
* Use Palabos vti writer;
* Use paraview to analyze the results;
* Use the voxelizer;
* Generate TriangleSets;
* Read external .stl files;
* Use sponge zones;
* Use off-lattice boundary conditions implemented in Palabos;
* Compute the drag coefficient.

**What it does not cover:**

* Moving objects
* Refilling algorithms

**Bibliography**

Bouzidi, M., Firdaouss, M., & Lallemand, P. (2001). Momentum transfer of a Boltzmann-lattice fluid with boundaries. Physics of Fluids, 13(11), 3452–3459. https://doi.org/10.1063/1.1399290

Yu, D., Mei, R., & Shyy, W. (2003). A Unified Boundary Treatment in Lattice Boltzmann Method. 41st Aerospace Sciences Meeting and Exhibit. https://doi.org/10.2514/6.2003-953

Ginzburg, I., Verhaeghe, F., & d’Humières, D. (2008). Two-Relaxation-Time Lattice Boltzmann Scheme: About Parametrization, Velocity, Pressure and Mixed Boundary Conditions. Commun. Comput. Phys., 3(2), 427–478.

Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., & Latt, J. (2021). Enhanced single-node lattice Boltzmann boundary condition for fluid flows. Physical Review E, 103(5), 053308. https://doi.org/10.1103/PhysRevE.103.053308

Marson, F. (2022). Directional lattice Boltzmann boundary conditions. https://doi.org/10.13097/archive-ouverte/unige:160770


# Guide

## 1. Install CMake, ParaView, and update your `c++` compiler

To compile this example you need CMake installed in your system. To check the installation type on your
terminal: `cmake --version`. If it returns an error, you need to proceed with the installation. Try with the following
command:

* Archlinux based (Manjaro): `sudo pacman -S cmake`
* Debian based (Ubuntu): `sudo apt install cmake`
* redhat based (OpenSuSE, Fedora, Centos): `yum install cmake`

for other systems, please refer to the [official installation instructions](https://cmake.org/install/).

### Update the c++ compiler

To compile this exercise, you need an updated compiler compatible with the c++17 standard.

* if you are using an Arch-based distribution, you are updated by definition.
* if you are using ubuntu, check this link https://linuxize.com/post/how-to-install-gcc-compiler-on-ubuntu-18-04/ and read the section "Installing Multiple GCC Versions".

### [optional, takes time, do it later] Install meshlab

Meshlab is an open-source program to manipulate meshes. It is particularly useful to scale and rotate `.stl` files before
importing them into Palabos. The easiest way is to install it from snap: `snap install meshlab`. Otherwise, you can google
for your distribution-dependent instructions.

## 2. Test the compilation

If Palabos is on the directory of the exercises, the provided CMakeLists.txt should already work. In all the other
cases, you can check the following lines of the CMakeLists.txt to modify the palabos location:

```cmake
# # # ADD PALABOS PATH:
# # CASE 1: PALABOS_ROOT is on your path
# file(TO_CMAKE_PATH $ENV{PALABOS_ROOT} PALABOS_ROOT)

# # CASE 2: specify an ABSOLUTE DIRECTORY PATH
#set(PALABOS_ROOT /home/mars/MEGA2/projects/palabos)

# # CASE 3: specify a relative path to the build folder (eg if you are in the showCases folder)
set(PALABOS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
```

Open a terminal in the current directory and type the following commands to compile the code with CMake:

```bash
cd build
cmake .. && make -j 2
cd ..
```

in the case of a Linux-based system. CMake runs in LIgeneric in release mode, if you want to change the build mode to
debug, you can use the CMake with the option `-DCMAKE_BUILD_TYPE=Debug`.

## 3. Compile and run the exercise after choosing the simulation parameters

Check now the file `curvedBoundaries_3d.cpp`. First, choose your numerical resolution in:

```c++
    IncomprFlowParam<Real> parameters =
        lu.initLrefluNodim(000000/*resolution*/, &dimless, 00000 /*u_lb*/,
                           false /*tau, optional*/)
            .printParameters()
            .getIncomprFlowParam();
```

then, compile it again, and run the simulation in parallel using mpirun (if installed in your system)

```bash
mpirun -np 2 curvedBoundaries_3d
```
in case of problems with mpi, try to run the serial version
```bash
./curvedBoundaries_3d
```

## 4. Drag coefficient

Now it is time to compute the drag coefficient, I can give you three hints:

1. Check the available methods of `OffLatticeBoundaryCondition3D<Real, Descriptor, Array<Real, 3>>`
2. Use the definition of [Drag coefficient](https://en.wikipedia.org/wiki/Drag_coefficient);
3. Use `pcout` to write output files in parallel simulations.

## 5. Warm-up

### Change boundary condition

At point 7 of the function `setupLattice3d()` you will find that a function is called to set up the off-lattice boundary
condition. Try to figure out how to modify it to use a different boundary condition.

### Use sponge zones

At point 5. of `setupLattice3d()` you will find a call to a function to set up a sponge zone. Check inside this helper
function to understand how it works. Then, try to set up sponge zones and to see their effects on the simulation. For
instance, compare the time evolution of the drag coefficient with and without its implementation.

### Add a second sphere and transform them into ellipsoids

Using either Guo (GZS) or Bouzidi (BFL), add a second sphere behind the first one. Then, try to transform them into
ellipsoids. Hint: check file `shapes.h`.

## 6. Have fun

### Use an external .stl mesh file

Inside this folder, you will find a file called `curvedBoundaries_2021.stl`. I downloaded it from https://www.thingiverse.com/,
then I have scaled, rotated, and simplified it using `meshlab`. Try to read it using the available method in `shapes.h`,
then try to run simulations with different parameters.

### [optional, it takes time, do it later] Use your own .stl file

Look for a stl file on https://www.thingiverse.com/. Then using meshlab:

1. center on the origin of the axes using the translate filter
2. rotate it using the rotation filter
3. reduce the amount triangles using a simplifying filter.

### From planar Poiseuille to periodic flow + neumann external boundary conditions:
Go to the following function:
```c++
template <typename Real, template <typename U> class Descriptor>
auto inject_on_lattice_bc(MultiBlockLattice3D<Real, Descriptor>* target_lattice,
                          IncomprFlowParam<Real> const& parameters)
```

you will see that initial and boundary conditions are imposed with the help of "functors":

* ``PoiseuilleVelocity<Real>(parameters)``;
* ``PoiseuilleDensity<Real,Descriptor>(parameters)``;
* ``PoiseuilleVelocityAndDensity<Real, Descriptor>(parameters)``.

Look for equivalent function objects that allow imposing a constant velocity and density profile.

Use periodic boundaries on the z-normal boundary faces of the domain, Neumann boundary conditions on the y-normal, and use a uniform pressure and velocity initial condition.

## 7. Implement HW

Try to implement the HW boundary condition starting from the BFL (Bouzidi) implementation. To begin with do:

```bash
cd src
cp path/to/palabos/src/offLattice/bouzidiOffLatticeModel3D.h HWLatticeModel3D.h
cp path/to/palabos/src/offLattice/bouzidiOffLatticeModel3D.hh HWLatticeModel3D.hh
```

As firs step, rename all the "bouzidiOffLattice" naming to the respective "HWLattice". For example, from:

```c++
#ifndef BOUZIDI_OFF_LATTICE_MODEL_3D_H
#define BOUZIDI_OFF_LATTICE_MODEL_3D_H
```

to

```c++
#ifndef HW_LATTICE_MODEL_3D_H
#define HW_LATTICE_MODEL_3D_H
```

Also, search for all the occurrences of `BouzidiOffLatticeModel3D` in the two files and replace them with `HWLatticeModel3D`. Do
the analogous thing also in the `.hh` file. Then comment out the following `#include` directives in the `.hh` file:

```c++
#include "offLattice/bouzidiOffLatticeModel3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include <algorithm>
#include <vector>
#include <cmath>
```

and replace them with

```c++
#include "palabos3D.h"
#include "palabos3D.hh"  // explicit inclusion because it contains templates
```

Finally, navigate to `setup.h` and include the following preprocessor directives:

```c++
#include "HWLatticeModel3D.h"
#include "HWLatticeModel3D.hh"
```

Now you are ready to implement the HW. Open `HWLatticeModel3D.hh` and navigate to

```c++
template<typename T, template<typename U> class Descriptor>
void HWLatticeModel3D<T,Descriptor>::cellCompletion ()
```

and try to modify the algorithm to use a simple half-way bounce-back. Finally, test it!


### Bug fix:
The following lines have been removed (they are not needed in this example, they slow the computations, but they should not change the results):
```c++
lattice.executeInternalProcessors();  // i.e. boundary conditions in
                                              // this case
lattice.incrementTime();
```
credits: Julius Weinmiller (Ulm University)



