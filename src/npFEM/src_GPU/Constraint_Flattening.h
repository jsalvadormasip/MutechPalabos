///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
#ifndef CONSTRAINT_FLATTENING_H
#define CONSTRAINT_FLATTENING_H
///////////////////////////////////////////////////////////////////////////////
#include "npFEM/src_GPU/common.h"
///////////////////////////////////////////////////////////////////////////////
namespace plb {
namespace npfem {
///////////////////////////////////////////////////////////////////////////////
struct Constraint_flat 
{
  int ConstraintType_ = -1;
  int idO_ = -1;
  ShapeOpScalar rangeMin_ = 1;
  ShapeOpScalar rangeMax_ = 1;
  ShapeOpScalar Scalar1_ = 1;
  ShapeOpScalar weight_ = 100;
  ShapeOpScalar E_nonePD_ = 0;
  int    *idI_;
  ShapeOpScalar *vectorx_;
  // Matrices are stored in ColMajor order which is also the default of Eigen
  ShapeOpScalar *matrix22_;
  ShapeOpScalar *matrix33_;


  /*
  int    idI_[4];
  ShapeOpScalar vectorx_[4] = {0,0,0,0};
  // Matrices are stored in ColMajor order which is also the default of Eigen
  ShapeOpScalar matrix22_[4] = { 0,0,0,0 };
  ShapeOpScalar matrix33_[9] = { 0,0,0, 0,0,0 ,0,0,0 };
  */
  
  Constraint_flat() 
  {
    ConstraintType_ = -1;
    idO_            = -1;
    rangeMin_       =  1.;
    rangeMax_       =  1.;
    Scalar1_        =  1.;
    weight_         =  100.;
    E_nonePD_       =  0.;
    idI_            = new int    [4];
    vectorx_        = new ShapeOpScalar [4];
    matrix22_       = new ShapeOpScalar [4];
    matrix33_       = new ShapeOpScalar [9];
    
    idI_[0] = -1     ; idI_[1] = -1     ; idI_[2] = -1     ; idI_[3] = -1     ;
    vectorx_[0]  = 0.; vectorx_[1]  = 0.; vectorx_[2]  = 0.; vectorx_[3]  = 0.;
    matrix22_[0] = 0.; matrix22_[1] = 0.; matrix22_[2] = 0.; matrix22_[3] = 0.;
    matrix33_[0] = 0.; matrix33_[1] = 0.; matrix33_[2] = 0.; matrix33_[3] = 0.; 
    matrix33_[4] = 0.; matrix33_[5] = 0.; matrix33_[6] = 0.; matrix33_[7] = 0.; 
    matrix33_[8] = 0.;
  }

  ~Constraint_flat() 
  {
    delete idI_;
    delete vectorx_;
    delete matrix22_;
    delete matrix33_;
  }
  
};
///////////////////////////////////////////////////////////////////////////////
} // namespace npfem
} // namespace plb
///////////////////////////////////////////////////////////////////////////////
#endif // CONSTRAINT_FLATTENING_H
///////////////////////////////////////////////////////////////////////////////