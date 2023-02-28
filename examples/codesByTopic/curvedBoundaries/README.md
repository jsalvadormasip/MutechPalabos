---
title: "Link-wise (directional) boundary conditions for LBM"
author: Francesco Marson, Department of Computer Science, University of Geneva
date: Palabos Summer School 2021
geometry: margin=2cm
output: pdf_document
---

The main goal of this example is to show how you can try and test different curved boundary conditions.
In particular, it shows how to use the the genericELI class that allows to define any LI+ methods.

LI+ methods includes two types of boundaries (Ginzburg et al. 2023):
* LI (e.g. Bouzidi et al. 2001, or Yu et al. 2003, or CLI Ginzburg et al. 2008)
* ELI: a class of local methods introduced in (Marson et al. 2021) and extended in (Marson 2022 and Ginzburg et al. 2023)

As beneficial collateral effects, you may familiarize yourself with the following topics

* Compile an example with CMake;
* Use Palabos vti writer;
* Use the voxelizer;
* Generate TriangleSets;
* Read external .stl files;
* Compute the drag coefficient.

**What it does not cover:**

* Moving objects
* Refilling algorithms

**Bibliography**

Filippova, O., & Hänel, D. (1997). Lattice-Boltzmann simulation of gas-particle flow in filters. Computers & Fluids, 26(7), 697–712. https://doi.org/10.1016/S0045-7930(97)00009-1

Bouzidi, M., Firdaouss, M., & Lallemand, P. (2001). Momentum transfer of a Boltzmann-lattice fluid with boundaries. Physics of Fluids, 13(11), 3452–3459. https://doi.org/10.1063/1.1399290

Guo, Z., Zheng, C., & Shi, B. (2002). An extrapolation method for boundary conditions in lattice Boltzmann method. Physics of Fluids, 14(6), 2007–2010. https://doi.org/10.1063/1.1471914


Yu, D., Mei, R., & Shyy, W. (2003). A Unified Boundary Treatment in Lattice Boltzmann Method. 41st Aerospace Sciences Meeting and Exhibit. https://doi.org/10.2514/6.2003-953

Ginzburg, I., Verhaeghe, F., & d’Humières, D. (2008). Two-Relaxation-Time Lattice Boltzmann Scheme: About Parametrization, Velocity, Pressure and Mixed Boundary Conditions. Commun. Comput. Phys., 3(2), 427–478.

Marson, F., Thorimbert, Y., Chopard, B., Ginzburg, I., & Latt, J. (2021). Enhanced single-node lattice Boltzmann boundary condition for fluid flows. Physical Review E, 103(5), 053308. https://doi.org/10.1103/PhysRevE.103.053308

Marson, F. (2022). Directional lattice Boltzmann boundary conditions. https://doi.org/10.13097/archive-ouverte/unige:160770

Ginzburg, I., Silva, G., Marson, F., Chopard, B., & Latt, J. (2023). Unified directional parabolic-accurate lattice Boltzmann boundary schemes for grid-rotated narrow gaps and curved walls in creeping and inertial fluid flows. Physical Review E, 107(2), 025303. https://doi.org/10.1103/PhysRevE.107.025303



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

## 3. Run 
```bash
mpirun -np 2 curvedBoundaries output_folder_name #METHOD#
```
where the possible values of `#METHOD#` are:
* `HW` (Ginzburg 1994)
* `BFL`(Bouzidi et al., 2001)
* `FH` (Filippova and Hanel, 1997)
* `MLS` (Mei, Luo, Shyy, 1999)
* `ELIULC`(Marson et al. 2021)
* `ELIUL` (Marson et al. 2021)
* `ELIULK1` (Marson et al. 2021)
* `ELIULK3` (Marson 2022, Ginzburg et al. 2023)
* `ELIULK4` (Marson 2022, Ginzburg et al. 2023)
* `LIgeneric` (Marson 2022, Ginzburg et al. 2023)

When `LIgeneric` is selected with for example 
```bash
mpirun -np 2 curvedBoundaries tmp_out_folder LIgeneric
```
the program will read the xml file `LIplusGeneric.xml` containing the following lines:
```xml
<!--admitted symbols: tauPlus, tauMinus, q, up, down -->
<generic_eli>
<alphaPlus>-1</alphaPlus>
<alphaMin>+1</alphaMin>
<beta>0</beta> <!--    if beta = 0 => ELI \in LI+, if beta \neq 0 => LI+-->
<kPlus>q - tauPlus</kPlus>
<kMin>0</kMin>
</generic_eli>
```
In this file you can change the values of the boundary condition parameters to match one of the infinite number of possible LI+ boundary conditions [but mind stability, see (Ginzburg et al., 2023)]. You can for example recover all the other pre-defined schemes in the list, or also other, inserting the values for alpha-, alpha+, beta, K+ and K- you can find, for example, in tables 5.3 pag 97 and 5.4 at page 102 of (Marson 2022). The above mentioned tables don't contain the values for K+ and K- which can be found 
* for LI in equaitons (5.8a, 5.13, 5.14)
* for ELI in equations (5.19) page 102.
you can find a more in-depth discussion of those methods in (Ginzburg et al. 2023).

in case of problems with mpi, try to run the serial version
```bash
./curvedBoundaries output_folder_name #METHOD#
```
Note: also the method of (Guo et al. 2002) is available in palabos and can be used in ths ecample updating the switch case in setup.h using `GuoOffLatticeModel3D<T, DESCRIPTOR>`
