cmake_minimum_required(VERSION 3.0)

find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    message("ccache used.")
endif() 

project(palabos)

enable_language(CXX)

###############################################################################
# Palabos Library
###############################################################################

file(GLOB_RECURSE PALABOS_SRC "src/*.cpp")
file(GLOB_RECURSE EXT_SRC "externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})
# PUBLIC directories are propagated to projects that links to palabos library
target_include_directories(palabos
    PUBLIC 
        "src"
        "externalLibraries"
        "externalLibraries/Eigen3"
)

# Add character 'd' to target compiled with DEBUG configuration
set_target_properties(palabos PROPERTIES DEBUG_POSTFIX "d")

###############################################################################

# This compiles Palabos library and all the examples
# It is used mainly for Continuous Integration (CI)
option(CI "Compile Everything (Continuous Integration)" OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
endif()
message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

# Compiler flags
# Append flags: set(CMAKE_XXX_FLAGS "${CMAKE_XXX_FLAGS} ...")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message("GCC.")
    target_compile_options(palabos PRIVATE 
                    -std=c++11 -Wall -Wnon-virtual-dtor
                    $<$<CONFIG:Release>:-O3 -DNDEBUG>
                    $<$<CONFIG:Debug>:-g -DPLB_DEBUG -O0>
    )
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Clang.")
    target_compile_options(palabos PRIVATE 
                    -std=c++11 -Wall -Wnon-virtual-dtor
                    $<$<CONFIG:Release>:-O3 -DNDEBUG>
                    $<$<CONFIG:Debug>:-g -DPLB_DEBUG -O0>
    )
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
    target_compile_options(palabos PRIVATE 
                    $<$<CONFIG:Release>:/Ox /Ot /GS- /GL /DNDEBUG>
                    $<$<CONFIG:Debug>:/DPLB_DEBUG>
    )
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message( FATAL_ERROR "CXX compiler not recognized. CMake will quit." )
endif()

option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message("Enabling MPI")
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        #set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        target_include_directories(palabos PUBLIC ${MPI_CXX_INCLUDE_PATH})
        # This function append to previous options
        target_compile_options(palabos PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${MPI_CXX_COMPILE_FLAGS}>)
        target_compile_definitions(palabos PUBLIC -DPLB_MPI_PARALLEL)
        # palabos does not need to link MPI to compile, however any code that consumes it must link MPI libraries
        target_link_libraries(palabos INTERFACE ${MPI_CXX_LIBRARIES})
    endif()
endif()

if(WIN32)
    option(ENABLE_POSIX "Enable POSIX" OFF)
else()
    option(ENABLE_POSIX "Enable POSIX" ON)
endif()

if(ENABLE_POSIX)
    message("Enabling POSIX")
    target_compile_definitions(palabos PRIVATE -DPLB_USE_POSIX)
endif()

if(APPLE)
    target_compile_definitions(palabos PRIVATE -DPLB_MAC_OS_X)
endif()

if(WIN32 OR CYGWIN)
    target_compile_definitions(-DPLB_WINDOWS)
endif()
###############################################################################

###############################################################################

# To link palabos in directories other than this
cmake_policy(SET CMP0079 NEW)

###############################################################################
# examples directory
###############################################################################
if(CI)
    # benchmarks
    add_subdirectory(examples/benchmarks/cavity3d)
    # codesByTopic
    add_subdirectory(examples/codesByTopic/asinariDynamics)
    add_subdirectory(examples/codesByTopic/bounceBack)
    add_subdirectory(examples/codesByTopic/boundaryCondition)
    add_subdirectory(examples/codesByTopic/complexNumbers)
    add_subdirectory(examples/codesByTopic/coProcessor)
    add_subdirectory(examples/codesByTopic/copying)
    add_subdirectory(examples/codesByTopic/couplings)
    add_subdirectory(examples/codesByTopic/cutTriangularSurface)
    add_subdirectory(examples/codesByTopic/d3q13)
    add_subdirectory(examples/codesByTopic/dataAnalysis)
    add_subdirectory(examples/codesByTopic/dotList)
    add_subdirectory(examples/codesByTopic/dynamicsInfo)
    add_subdirectory(examples/codesByTopic/externalMacroscopicVariables)
    add_subdirectory(examples/codesByTopic/io)
    add_subdirectory(examples/codesByTopic/linearAlgebra)
    add_subdirectory(examples/codesByTopic/marchingCube)
    add_subdirectory(examples/codesByTopic/multiBlock)
    add_subdirectory(examples/codesByTopic/navierStokesModels)
    add_subdirectory(examples/codesByTopic/particlesInTube)
    add_subdirectory(examples/codesByTopic/reparallelize)
    add_subdirectory(examples/codesByTopic/scalarField)
    add_subdirectory(examples/codesByTopic/shanChenMultiPhase)
    add_subdirectory(examples/codesByTopic/smagorinskyModel)
    add_subdirectory(examples/codesByTopic/transientStatistics)
    add_subdirectory(examples/codesByTopic/userInput)
    add_subdirectory(examples/codesByTopic/aneurysm_xdmf_hdf5_output)
    add_subdirectory(examples/codesByTopic/TRTdynamics)
    # showCases
    add_subdirectory(examples/showCases/aneurysm)
    add_subdirectory(examples/showCases/bloodFlowDefoBodies)
    add_subdirectory(examples/showCases/boussinesqThermal2d)
    add_subdirectory(examples/showCases/boussinesqThermal3d)
    add_subdirectory(examples/showCases/carreauPoiseuille)
    add_subdirectory(examples/showCases/cavity2d)
    add_subdirectory(examples/showCases/cavity3d)
    add_subdirectory(examples/showCases/collidingBubbles3d)
    add_subdirectory(examples/showCases/cylinder2d)
    add_subdirectory(examples/showCases/cylinder3d)
    add_subdirectory(examples/showCases/externalFlowAroundObstacle)
    add_subdirectory(examples/showCases/generalExternalFlow)
    add_subdirectory(examples/showCases/gridRefinement2d)
    add_subdirectory(examples/showCases/gridRefinement3d)
    add_subdirectory(examples/showCases/heatedObjectInChannel)
    add_subdirectory(examples/showCases/movingWall)
    add_subdirectory(examples/showCases/multiComponent2d)
    add_subdirectory(examples/showCases/multiComponent3d)
    add_subdirectory(examples/showCases/particlesInCone)
    add_subdirectory(examples/showCases/poiseuille)
    add_subdirectory(examples/showCases/rectangularChannel3d)
    add_subdirectory(examples/showCases/vofMultiPhase)
    add_subdirectory(examples/showCases/womersley)
    # tutorial
    add_subdirectory(examples/tutorial/permeability)
    add_subdirectory(examples/tutorial/tutorial_1)
    add_subdirectory(examples/tutorial/tutorial_2)

    ###############################################################################
    # utility directory
    ###############################################################################
    # stl
    add_subdirectory(utility/stl/cutSTL)
    add_subdirectory(utility/stl/mergeSTL)
    add_subdirectory(utility/stl/refineByAreaSTL)
    add_subdirectory(utility/stl/refineSTL)
    add_subdirectory(utility/stl/reverseOrientationSTL)
    add_subdirectory(utility/stl/revoxelize)
    add_subdirectory(utility/stl/rotateAtOriginSTL)
    add_subdirectory(utility/stl/scaleDirectionallySTLraw)
    add_subdirectory(utility/stl/scaleSTL)
    add_subdirectory(utility/stl/scaleSTLraw)
    add_subdirectory(utility/stl/shiftSTLraw)
    add_subdirectory(utility/stl/stlCheck)
    add_subdirectory(utility/stl/stlInfo)
    add_subdirectory(utility/stl/stlProjectedArea)
    add_subdirectory(utility/stl/stlShadowArea)
    add_subdirectory(utility/stl/stlToHtml)
    add_subdirectory(utility/stl/stlToX3D)
    add_subdirectory(utility/stl/toAsciiSTL)
    add_subdirectory(utility/stl/toBinarySTL)
    # vtk
    add_subdirectory(utility/vtk/plb2vtk2D)
    add_subdirectory(utility/vtk/plb2vtk3D)
endif()