///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifndef LSSOLVER_H
#define LSSOLVER_H
///////////////////////////////////////////////////////////////////////////////
#include "Types.h"
//#include <unsupported/Eigen/src/IterativeSolvers/MINRES.h>
///////////////////////////////////////////////////////////////////////////////
/* This file contains all the linear system solvers of the ShapeOp library. */
///////////////////////////////////////////////////////////////////////////////
namespace ShapeOp {
///////////////////////////////////////////////////////////////////////////////
/* 
  Base class of any sparse linear system solver. This class defines the
  main functionalities of the ShapeOp sparse linear system solvers (Ax = b).
*/
class SHAPEOP_API LSSolver {
public:
    virtual ~LSSolver(){};
    /* Initialize the linear system solver using the sparse matrix A. */
    virtual void initialize(const SparseMatrix& A, unsigned int iteration = 1)
        = 0;
    /* Solve the linear system Ax = b. */
    virtual VectorX solve(const VectorX& b, const VectorX& x0) const = 0;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const = 0;
};
///////////////////////////////////////////////////////////////////////////////
/* 
Sparse linear system solver based on Cholesky. This class implements
a sparse linear system solver based on the Cholesky LDL^T algorithm from
Eigen.
*/
class SHAPEOP_API SimplicialLDLTSolver : public LSSolver {
public:
    virtual ~SimplicialLDLTSolver(){};
    /* Prefactorize the sparse matrix (A = LDL^T). */
    virtual void initialize(
        const SparseMatrix& A, unsigned int iteration) override final;
    /* Solve the linear system by applying twice backsubstitution. */
    virtual VectorX solve(
        const VectorX& b, const VectorX& x0) const override final;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const override final;
private:
    Eigen::SimplicialLDLT<SparseMatrix> solver_;
};
///////////////////////////////////////////////////////////////////////////////
/*
Sparse linear system solver based on CG. This class implements a
sparse linear system solver based on the CG algorithm from Eigen.
*/
class SHAPEOP_API CGSolver : public LSSolver {
public:
    virtual ~CGSolver(){};
    /* Initialize PCG. */
    virtual void initialize(
        const SparseMatrix& A, unsigned int iteration) override final;
    /* Solve the linear system by applying CG. */
    virtual VectorX solve(
        const VectorX& b, const VectorX& x0) const override final;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const override final;
private:
    Eigen::ConjugateGradient<SparseMatrix, Eigen::Lower,
        Eigen::IncompleteLUT<Scalar>>
        solver_;
};
///////////////////////////////////////////////////////////////////////////////
/*
Sparse linear system solver based on MINRES. This class implements a
sparse linear system solver based on the MINRES algorithm from Eigen.

/*
class SHAPEOP_API MINRESSolver : public LSSolver {
public:
    virtual ~MINRESSolver(){};
    // Initialize MINRES.
    virtual void initialize(
    const SparseMatrix& A, unsigned int iteration) override final;
    //Solve the linear system by applying MINRES.
    virtual VectorX solve(
    const VectorX& b, const VectorX& x0) const override final;
    // Reports whether previous computation was successful.
    virtual Eigen::ComputationInfo info() const override final;
private:
    Eigen::MINRES<SparseMatrix, Eigen::Lower, Eigen::IncompleteLUT<Scalar>>
    solver_;
};
*/
///////////////////////////////////////////////////////////////////////////////
/* Sparse linear system solver based on successive over-relaxation (SOR). */
class SHAPEOP_API SORSolver : public LSSolver {
public:
    SORSolver(ShapeOp::Scalar relaxation = 1.6);
    virtual ~SORSolver(){};
    /* Initialize SOR. */
    virtual void initialize(
        const SparseMatrix& A, unsigned int iteration) override final;
    /* Solve the linear system by applying SOR. */
    virtual VectorX solve(
        const VectorX& b, const VectorX& x0) const override final;
    /* Reports whether previous computation was successful. */
    virtual Eigen::ComputationInfo info() const override final;
private:
    SparseMatrixT<Eigen::RowMajor> A_;
    ShapeOp::Scalar relaxation_;
    unsigned int iteration_;
};
///////////////////////////////////////////////////////////////////////////////
} // namespace ShapeOp
///////////////////////////////////////////////////////////////////////////////
#ifdef SHAPEOP_HEADER_ONLY
#include "LSSolver.cpp"
#endif
///////////////////////////////////////////////////////////////////////////////
#endif // LSSOLVER_H
///////////////////////////////////////////////////////////////////////////////
