# Instructions for npFEM (ShapeOp) with Rhino3D-Grasshopper

This is a stand-alone version of the library, in the sense that it is completely detached from Palabos.

Whatever presented here is intended for Windows only because Rhino3D & Grasshopper (https://www.rhino3d.com/) are available only for Windows (and Mac recently). Also, we are using Rhino version 5 and Grasshopper (GH) installed separately (in latest Rhino, GH is integrated in it). The last choice is because we want to keep this code compatible with the ShapeOp project (https://www.shapeop.org/).

npFEM is a heavily modified version of ShapeOp.

For more information about ShapeOp, and how it can be integrated in Rhino-GH, you can consult the links below:

* https://www.shapeop.org/
* https://doi.org/10.1007/978-3-319-24208-8_42

Essentially, we are going to compile the library into a dynamic library (dll), and use Rhino-GH for simple tasks and visualization. By simple tasks, we mean one RBC/PLT only, as this is intended for building the materials and then feed them into the Palabos-npFEM library for more intricate simulations.

This part of the project is inspired by the project below:

* https://github.com/AndersDeleuran/ShapeOpGHPython

What we provide in the npFEM_RhinoGH folder adopts ideas and code from the project above, and it is tailored for blood cells.

Here, we simply modify ShapeOp (npFEM) and use Rhino-GH for preliminary steps.

## Compilation

Windows only (CMake & Visual Studio):

1. mkdir build
2. cd build
3. cmake ..
4. make -j

This results in a dynamic library inside the build folder.

* Follow the instructions at ShapeOp page on how to install the Rhino-GH environment.
* Copy the produced dll in the libraries folder of Grasshopper.
* Open Rhino-GH and the RBC.gh (Grasshopper file) provided in the npFEM_RhinoGH folder.

## Rhino-GH

* Go to npFEM_RhinoGH folder
* Open the RBC.3dm Rhino File (this will open Rhino3D)
* In Rhino's command-line type: Grasshopper (this will open GH)
* In GH's environment open RBC.gh file
* The rest are similar to what is described in this project: https://github.com/AndersDeleuran/ShapeOpGHPython
