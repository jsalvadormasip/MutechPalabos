///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Remark on RowMajor:
// I have tried to use RowMajor as well, but there is an assertion failure
// for column vectors. Essentially, column vectors are not allowed by Eigen
// to be RowMajor (which makes sense).
// There is a work-around by using the ternary operator in the Options, but
// for now we go with the Eigen default which is ColMajor
///////////////////////////////////////////////////////////////////////////////
#ifndef TYPES_H
#define TYPES_H
///////////////////////////////////////////////////////////////////////////////
#include "src_GPU/common.h"
// Stand-alone development
#ifdef NO_PALABOS
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include "../../shapeOp_GH/externalLibraries/eigen/unsupported/Eigen/KroneckerProduct"
#else
// Eigen used by Palabos
#include <Eigen3/Core>
#include <Eigen3/Dense>
#include <Eigen3/Sparse>
//#include <Eigen3/Eigenvalues>
//#include "../../shapeOp_GH/externalLibraries/eigen/unsupported/Eigen/KroneckerProduct"
#endif
// For Collisions
#ifndef NO_PALABOS
#include "nanoflann/nanoflann.hpp"
#else
#include "nanoflann.hpp"
#endif
///////////////////////////////////////////////////////////////////////////////
/**
This file redefines EIGEN types using the scalar type ::ShapeOpScalar defined in
Common.h.*/
///////////////////////////////////////////////////////////////////////////////
/** Defines Eigen Alignment type.*/
#ifdef SHAPEOP_DONT_ALIGN
#define SHAPEOP_ALIGNMENT Eigen::DontAlign
#else
#define SHAPEOP_ALIGNMENT Eigen::AutoAlign
#endif
///////////////////////////////////////////////////////////////////////////////
/** \brief Namespace of the ShapeOp library.*/
namespace ShapeOp {
typedef ShapeOpScalar Scalar; // A scalar type, double or float, as defined in
// ::ShapeOpScalar in Common.h.
// Dense
// Biwise OR |: In our case is 0 | 0 => 0 => ColMajor
template <int Rows, int Cols,
    int Options = (Eigen::ColMajor | SHAPEOP_ALIGNMENT)>
using MatrixT = Eigen::Matrix<Scalar, Rows, Cols,
    Options>; // A typedef of the dense matrix of Eigen.
typedef MatrixT<2, 1> Vector2; // A 2d column vector.
typedef MatrixT<2, 2> Matrix22; // A 2 by 2 matrix.
typedef MatrixT<2, 3> Matrix23; // A 2 by 3 matrix.
typedef MatrixT<3, 1> Vector3; // A 3d column vector.
typedef MatrixT<3, 2> Matrix32; // A 3 by 2 matrix.
typedef MatrixT<3, 3> Matrix33; // A 3 by 3 matrix.
typedef MatrixT<3, 4> Matrix34; // A 3 by 4 matrix.
typedef MatrixT<4, 1> Vector4; // A 4d column vector.
typedef MatrixT<4, 4> Matrix44; // A 4 by 4 matrix.
typedef MatrixT<3, Eigen::Dynamic> Matrix3X; // A 3 by n matrix.
typedef MatrixT<Eigen::Dynamic, 3> MatrixX3; // A n by 3 matrix.
typedef MatrixT<Eigen::Dynamic, 1> VectorX; // A nd column vector.
typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic> MatrixXX; // A n by m matrix.
typedef MatrixT<Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXXCuda; // A n by m matrix.
// Sparse
template <int Options = Eigen::ColMajor>
using SparseMatrixT = Eigen::SparseMatrix<Scalar,
    Options>; // A typedef of the sparse matrix of Eigen.
typedef SparseMatrixT<> SparseMatrix; // The default sparse matrix of Eigen.
typedef Eigen::Triplet<Scalar>
    Triplet; // A triplet, used in the sparse triplet representation for
// matrices.
///////////////////////////////////////////////////////////////////////////////
} // namespace ShapeOp
///////////////////////////////////////////////////////////////////////////////
#endif // TYPES_H
///////////////////////////////////////////////////////////////////////////////
