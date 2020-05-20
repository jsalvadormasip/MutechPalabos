///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifndef COMMON_H
#define COMMON_H
///////////////////////////////////////////////////////////////////////////////
/**
This file is used to define macros and typedefs used by the C++ code and the C
API.*/
///////////////////////////////////////////////////////////////////////////////
/** Defines the scalar type used by the ShapeOp solver (float or double).*/
#ifndef SHAPEOP_SCALAR
// Change both defs
// Depending on the decision the solver may be more unstable
// All the testing is done with double
#define SHAPEOP_SCALAR double // double or float
#define DOUBLE_SHAPEOP // {FLOAT or DOUBLE}_SHAPEOP
#endif

#define MAX_ATTEMPT 7
#define MAX_LINE_SEARCH 5
#define GAMMA 0.0001
#define GAMMA2 0.9
#define TOL 0.00001
#define CONV_WINDOWS 5
#define MEM_SIZE 5

#define BL 0

/** Defines the scalar type used by the ShapeOp solver (float or double).*/
typedef SHAPEOP_SCALAR ShapeOpScalar;
typedef double cuda_scalar;
///////////////////////////////////////////////////////////////////////////////
#define LU_GPU // Faster than the other 2
//#define CHOLESKY_GPU
//#define QR_GPU //By far slower
///////////////////////////////////////////////////////////////////////////////
/** Defines the API prefix for the current platform.*/
#if defined(_WIN32) || defined(_WIN64)
#pragma warning(disable : 4251) // Disable warnings about templates and std
                                // types exposed in the c++ interface.
#pragma warning(disable : 4244) // conversion from 'double' to
                                // 'ShapeOp::Scalar', possible loss of data
#pragma warning(disable : 4996) // dll window function
#ifdef SHAPEOP_EXPORT
#define SHAPEOP_API __declspec(dllexport)
#else
#define SHAPEOP_API __declspec(dllimport)
#endif
#else
#define SHAPEOP_API
#endif
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#define SHAPEOP_INLINE inline
#else
#define SHAPEOP_INLINE
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // COMMON_H
///////////////////////////////////////////////////////////////////////////////
