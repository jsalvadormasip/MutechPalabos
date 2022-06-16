#ifndef ADVECTION_DIFFUSION_FD_HH
#define ADVECTION_DIFFUSION_FD_HH

#include "WENO3_procedure.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/dynamics.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"

#include <cmath>

namespace plb {


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//3rd order Weighted Essentially Non-Oscillatory
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
WENO3<T>::WENO3(T d_, T eps_, bool neumann_, plint nx_, plint ny_, plint nz_)
: d(d_), eps(eps_), neumann(neumann_), nx(nx_), ny(ny_), nz(nz_)
{ }

template<typename T>
void WENO3<T>::processGenericBlocks (
        Box3D domain,  std::vector<AtomicBlock3D*> fields )
{
  PLB_PRECONDITION( fields.size()==6 );
  ScalarField3D<T>* phi_t    = dynamic_cast<ScalarField3D<T>*>(fields[0]);
  ScalarField3D<T>* phi_tp1  = dynamic_cast<ScalarField3D<T>*>(fields[1]);
  ScalarField3D<T>* result   = dynamic_cast<ScalarField3D<T>*>(fields[2]);
  TensorField3D<T,3>* uField = dynamic_cast<TensorField3D<T,3>*>(fields[3]);
  ScalarField3D<T>* Q        = dynamic_cast<ScalarField3D<T>*>(fields[4]);
  ScalarField3D<T>* v_sedimentation       = dynamic_cast<ScalarField3D<T>*>(fields[5]);



  Dot3D absoluteOffset = phi_tp1->getLocation();

  Dot3D ofs1 = computeRelativeDisplacement(*phi_t, *phi_tp1);
  Dot3D ofs2 = computeRelativeDisplacement(*phi_t, *result);
  Dot3D ofs3 = computeRelativeDisplacement(*phi_t, *uField);
  Dot3D ofs4 = computeRelativeDisplacement(*phi_t, *Q);
  Dot3D ofs5 = computeRelativeDisplacement(*phi_t, *v_sedimentation);

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {

                  plint absoluteX = absoluteOffset.x + iX + ofs1.x;
                  plint absoluteY = absoluteOffset.y + iY + ofs1.y;
                  plint absoluteZ = absoluteOffset.z + iZ + ofs1.z;



T phiC = phi_tp1->get(iX  +ofs1.x, iY  +ofs1.y, iZ  +ofs1.z);
T phiE1 = phi_tp1->get(iX+1+ofs1.x, iY  +ofs1.y, iZ  +ofs1.z);
T phiW1 = phi_tp1->get(iX-1+ofs1.x, iY  +ofs1.y, iZ  +ofs1.z);
T phiN1 = phi_tp1->get(iX  +ofs1.x, iY+1+ofs1.y, iZ  +ofs1.z);
T phiS1 = phi_tp1->get(iX  +ofs1.x, iY-1+ofs1.y, iZ  +ofs1.z);
T phiT1 = phi_tp1->get(iX  +ofs1.x, iY  +ofs1.y, iZ+1+ofs1.z);
T phiB1 = phi_tp1->get(iX  +ofs1.x, iY  +ofs1.y, iZ-1+ofs1.z);

T phiE2 = phi_tp1->get(iX+2+ofs1.x, iY  +ofs1.y, iZ  +ofs1.z);
T phiW2 = phi_tp1->get(iX-2+ofs1.x, iY  +ofs1.y, iZ  +ofs1.z);
T phiN2 = phi_tp1->get(iX  +ofs1.x, iY+2+ofs1.y, iZ  +ofs1.z);
T phiS2 = phi_tp1->get(iX  +ofs1.x, iY-2+ofs1.y, iZ  +ofs1.z);
T phiT2 = phi_tp1->get(iX  +ofs1.x, iY  +ofs1.y, iZ+2+ofs1.z);
T phiB2 = phi_tp1->get(iX  +ofs1.x, iY  +ofs1.y, iZ-2+ofs1.z);



Array<T,3> const& u = uField->get(iX+ofs3.x, iY+ofs3.y, iZ+ofs3.z);
T v_sedC = v_sedimentation->get(iX+ofs5.x, iY+ofs5.y, iZ+ofs5.z);
T v_sedT1 = v_sedimentation->get(iX+ofs5.x, iY+ofs5.y, iZ+1+ofs5.z);
T v_sedB1 = v_sedimentation->get(iX+ofs5.x, iY+ofs5.y, iZ-1+ofs5.z);


Array<T, 3> adv;

Array<T, 3> fp_p12;
Array<T, 3> fp_n12;

Array<T, 3> fp_p12_1;
Array<T, 3> fp_p12_2;
Array<T, 3> fp_n12_1;
Array<T, 3> fp_n12_2;

Array<T, 3> bp_p12_1;
Array<T, 3> bp_p12_2;
Array<T, 3> bp_n12_1;
Array<T, 3> bp_n12_2;

Array<T, 3> alpha_p_p12_1;
Array<T, 3> alpha_p_p12_2;
Array<T, 3> alpha_p_n12_1;
Array<T, 3> alpha_p_n12_2;

Array<T, 3> w1p_p12;
Array<T, 3> w2p_p12;
Array<T, 3> w1p_n12;
Array<T, 3> w2p_n12;


  if (util::greaterThan(u[0], (T) 0)) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_p12x
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  fp_p12_1[0] = (phiC + phiE1)/(T)2;
  fp_p12_2[0] = -(phiW1 - (T)3 * phiC)/(T)2;

  bp_p12_1[0] = (phiE1 - phiC)*(phiE1 - phiC);
  bp_p12_2[0] = (phiC - phiW1)*(phiC - phiW1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_n12x
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  fp_n12_1[0] = (phiW1 + phiC)/(T)2;
  fp_n12_2[0] = -(phiW2 - (T)3*phiW1)/(T)2;

  bp_n12_1[0] = (phiC - phiW1)*(phiC - phiW1);
  bp_n12_2[0] = (phiW1 - phiW2)*(phiW1 - phiW2);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


} else {

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //  fp_p12x
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    fp_p12_1[0] = -(phiE2 - (T)3*phiE1)/(T)2;
    fp_p12_2[0] = (phiE1 + phiC)/(T)2;

    bp_p12_1[0] = (phiE1 - phiE2)*(phiE1 - phiE2);
    bp_p12_2[0] = (phiC - phiE1)*(phiC - phiE1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_n12x
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  fp_n12_1[0] = -(phiE1 - (T)3*phiC)/(T)2;
  fp_n12_2[0] = (phiC + phiW1)/(T)2;

  bp_n12_1[0] = (phiC - phiE1)*(phiC - phiE1);
  bp_n12_2[0] = (phiW1 - phiC)*(phiW1 - phiC);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // adv[0] = fp_p12x - fp_n12x;


}

if (util::greaterThan(u[1], (T) 0)) {
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_p12y
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  fp_p12_1[1] = (phiC + phiN1)/(T)2;
  fp_p12_2[1] = -(phiS1 - (T)3 * phiC)/(T)2;


  bp_p12_1[1] = (phiN1 - phiC)*(phiN1 - phiC);
  bp_p12_2[1] = (phiC - phiS1)*(phiC - phiS1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_n12y
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fp_n12_1[1] = (phiS1 + phiC)/(T)2;
fp_n12_2[1] = -(phiS2 - (T)3*phiS1)/(T)2;


bp_n12_1[1] = (phiC - phiS1)*(phiC - phiS1);
bp_n12_2[1] = (phiS1 - phiS2)*(phiS1 - phiS2);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // adv[1] = fp_p12y - fp_n12y;
} else {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_p12y
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fp_p12_1[1] = -(phiN2 - (T)3*phiN1)/(T)2;
fp_p12_2[1] = (phiN1 + phiC)/(T)2;

bp_p12_1[1] = (phiN1 - phiN2)*(phiN1 - phiN2);
bp_p12_2[1] = (phiC - phiN1)*(phiC - phiN1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_n12y
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fp_n12_1[1] = -(phiN1 - (T)3*phiC)/(T)2;
fp_n12_2[1] = (phiC + phiS1)/(T)2;

bp_n12_1[1] = (phiC - phiN1)*(phiC - phiN1);
bp_n12_2[1] = (phiS1 - phiC)*(phiS1 - phiC);


}

if (util::greaterThan((u[2]+v_sedC), (T) 0)) {

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_p12z
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  fp_p12_1[2] = (phiC + phiT1)/(T)2;
  fp_p12_2[2] = -(phiB1 - (T)3 * phiC)/(T)2;


  bp_p12_1[2] = (phiT1 - phiC)*(phiT1 - phiC);
  bp_p12_2[2] = (phiC - phiB1)*(phiC - phiB1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_n12z
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fp_n12_1[2] = (phiB1 + phiC)/(T)2;
fp_n12_2[2] = -(phiB2 - (T)3*phiB1)/(T)2;


bp_n12_1[2] = (phiC - phiB1)*(phiC - phiB1);
bp_n12_2[2] = (phiB1 - phiB2)*(phiB1 - phiB2);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} else {


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_p12z
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fp_p12_1[2] = -(phiT2 - (T)3*phiT1)/(T)2;
fp_p12_2[2] = (phiT1 + phiC)/(T)2;

bp_p12_1[2] = (phiT1 - phiT2)*(phiT1 - phiT2);
bp_p12_2[2] = (phiC - phiT1)*(phiC - phiT1);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  fp_n12z
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fp_n12_1[2] = -(phiT1 - (T)3*phiC)/(T)2;
fp_n12_2[2] = (phiC + phiB1)/(T)2;

bp_n12_2[2] = (phiB1 - phiC)*(phiB1 - phiC);
bp_n12_1[2] = (phiC - phiT1)*(phiC - phiT1);

// /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}

  for(plint i = 0; i <= 2; i++){

  alpha_p_p12_1[i] = (T)2 / ((T)3 * (eps + bp_p12_1[i])*(eps + bp_p12_1[i]));
  alpha_p_p12_2[i] = (T)1 / ((T)3 * (eps + bp_p12_2[i])*(eps + bp_p12_2[i]));
  alpha_p_n12_1[i] = (T)2 / ((T)3 * (eps + bp_n12_1[i])*(eps + bp_n12_1[i]));
  alpha_p_n12_2[i] = (T)1 / ((T)3 * (eps + bp_n12_2[i])*(eps + bp_n12_2[i]));


  w1p_p12[i] = alpha_p_p12_1[i] / (alpha_p_p12_1[i] + alpha_p_p12_2[i]);
  w2p_p12[i] = alpha_p_p12_2[i] / (alpha_p_p12_1[i] + alpha_p_p12_2[i]);
  w1p_n12[i] = alpha_p_n12_1[i] / (alpha_p_n12_1[i] + alpha_p_n12_2[i]);
  w2p_n12[i] = alpha_p_n12_2[i] / (alpha_p_n12_1[i] + alpha_p_n12_2[i]);


  fp_p12[i] = w1p_p12[i] * fp_p12_1[i] + w2p_p12[i] * fp_p12_2[i];
  fp_n12[i] = w1p_n12[i] * fp_n12_1[i] + w2p_n12[i] * fp_n12_2[i];

}

  adv = fp_p12 - fp_n12;
  T diffX, diffY, diffZ;

  diffX = phiW1 + phiE1 - (T) 2 * phiC;
  diffY = phiS1 + phiN1 - (T) 2 * phiC;
  diffZ = phiT1 + phiB1 - (T) 2 * phiC;


if(neumann){

  if(absoluteX==0)
  {
    (util::greaterThan(u[0], (T) 0) ? (phiW1 = phiC) : (util::lessThan(   u[0], (T) 0) ? (phiC =  phiE1) : (phiW1 = phiE1)));
    adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));
    adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));
    adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));
    diffX = (T) 2 * phiE1 - (T) 2 * phiC;
  }

  if(absoluteX==nx-1)
  {
    (util::greaterThan(u[0], (T) 0) ? (phiC = phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 = phiC) : (phiE1 = phiW1)));
    adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));
    adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));
    adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));
    diffX = (T) 2 * phiW1 - (T) 2 * phiC;
  }

  if(absoluteY==0)
  {
    (util::greaterThan(u[1], (T) 0) ? (phiS1 = phiC) : (util::lessThan(   u[1], (T) 0) ? (phiC =  phiN1) : (phiS1 = phiN1)));
    adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));
    adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));
    adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));
        diffY = (T) 2 * phiN1 - (T) 2 * phiC;
  }

  if(absoluteY==ny-1)
  {
    (util::greaterThan(u[1], (T) 0) ? (phiC = phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 = phiC) : (phiN1 = phiS1)));
    adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));
    adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));
    adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));
    diffY = (T) 2 * phiS1 - (T) 2 * phiC;
  }

  if(absoluteZ==0)
  {
    (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiB1 = phiC) : (util::lessThan(   (u[2]+v_sedC), (T) 0) ? (phiC =  phiT1) : (phiB1 = phiT1)));
    adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));
    adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));
    adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));
    diffZ = (T) 2 * phiT1 - (T) 2 * phiC;
  }

  if(absoluteZ==nz-1)
  {
    (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC = phiB1) : (util::lessThan(   (u[2]+v_sedC), (T) 0) ? (phiT1 = phiC) : (phiT1 = phiB1)));
    adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));
    adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));
    adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));
    diffZ = (T) 2 * phiB1 - (T) 2 * phiC;
  }


  if(absoluteX==1) { adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));}
  if(absoluteX==1) { adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));}
  if(absoluteX==1) { adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));}


  if(absoluteX==nx-2) { adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));}
  if(absoluteX==nx-2) { adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));}
  if(absoluteX==nx-2) { adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));}


  if(absoluteY==1) { adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));}
  if(absoluteY==1) { adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));}
  if(absoluteY==1) { adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));}


  if(absoluteY==ny-2) { adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));}
  if(absoluteY==ny-2) { adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));}
  if(absoluteY==ny-2) { adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));}


  if(absoluteZ==1) { adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));}
  if(absoluteZ==1) { adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));}
  if(absoluteZ==1) { adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));}


  if(absoluteZ==nz-2) { adv[0] = (util::greaterThan(u[0], (T) 0) ? (phiC - phiW1) : (util::lessThan(   u[0], (T) 0) ? (phiE1 - phiC) : (T) 0.5 * (phiE1 - phiW1)));}
  if(absoluteZ==nz-2) { adv[1] = (util::greaterThan(u[1], (T) 0) ? (phiC - phiS1) : (util::lessThan(   u[1], (T) 0) ? (phiN1 - phiC) : (T) 0.5 * (phiN1 - phiS1)));}
  if(absoluteZ==nz-2) { adv[2] = (util::greaterThan((u[2]+v_sedC), (T) 0) ? (phiC - phiB1) : (util::lessThan((u[2]+v_sedC), (T) 0) ? (phiT1 - phiC) : (T) 0.5 * (phiT1 - phiB1)));}


  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//     DIFFUSION
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  T diff = d * (diffX + diffY + diffZ);


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


T advection = u[0] * adv[0] + u[1] * adv[1] + (u[2] + v_sedC) * adv[2];


T accu = phiC * (v_sedT1 - v_sedB1)/2;

result->get(iX+ofs2.x,iY+ofs2.y,iZ+ofs2.z) = - advection + Q->get(iX+ofs4.x,iY+ofs4.y,iZ+ofs4.z) - accu + diff;


                }
            }
        }
  }


template<typename T>
WENO3<T>* WENO3<T>::clone() const
{
  return new WENO3<T>(*this);
}

template<typename T>
void WENO3<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
  modified[0] = modif::nothing;           // phi_t
  modified[1] = modif::nothing;           // phi_tp1
  modified[2] = modif::staticVariables;   // result
  modified[3] = modif::nothing;           // u
  modified[4] = modif::nothing;           // Q
  modified[5] = modif::nothing;



}

/******************************3rd order Runge-Kutta for the time discretization***************************************************/
/* ******** RK3_Step1_functional3D ****************************************** */

template<typename T>
void RK3_Step1_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size() == 3 );
    ScalarField3D<T>& phi_n = *fields[0];
    ScalarField3D<T>& phi_n_adv = *fields[1];
    ScalarField3D<T>&  phi_1 = *fields[2];
    Dot3D offset_phi_n_adv  = computeRelativeDisplacement(phi_n, phi_n_adv);
    Dot3D offset_phi_1      = computeRelativeDisplacement(phi_n, phi_1);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                phi_1.get(iX+offset_phi_1.x,iY+offset_phi_1.y,iZ+offset_phi_1.z)
                    = phi_n.get(iX,iY,iZ) +
                      phi_n_adv.get(iX+offset_phi_n_adv.x,iY+offset_phi_n_adv.y,iZ+offset_phi_n_adv.z);
            }
        }
    }
}

template<typename T>
RK3_Step1_functional3D<T>* RK3_Step1_functional3D<T>::clone() const {
    return new RK3_Step1_functional3D<T>(*this);
}

template<typename T>
void RK3_Step1_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // phi_n
    modified[1] = modif::nothing; // phi_n_adv
    modified[2] = modif::staticVariables; // phi_1
}

template<typename T>
BlockDomain::DomainT RK3_Step1_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope; // Everything is local, no communication needed.
}


/* ******** RK3_Step2_functional3D ****************************************** */

template<typename T>
void RK3_Step2_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size() == 4 );
    ScalarField3D<T>& phi_n = *fields[0];
    ScalarField3D<T>& phi_1 = *fields[1];
    ScalarField3D<T>& phi_1_adv = *fields[2];
    ScalarField3D<T>& phi_2 = *fields[3];
    Dot3D offset_phi_1      = computeRelativeDisplacement(phi_n, phi_1);
    Dot3D offset_phi_1_adv  = computeRelativeDisplacement(phi_n, phi_1_adv);
    Dot3D offset_phi_2      = computeRelativeDisplacement(phi_n, phi_2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                phi_2.get(iX+offset_phi_2.x,iY+offset_phi_2.y,iZ+offset_phi_2.z)
                    = 3./4. * phi_n.get(iX,iY,iZ) +
                      1./4. * phi_1.get(iX+offset_phi_1.x,iY+offset_phi_1.y,iZ+offset_phi_1.z) +
                      1./4. * phi_1_adv.get(iX+offset_phi_1_adv.x,iY+offset_phi_1_adv.y,iZ+offset_phi_1_adv.z);
            }
        }
    }
}

template<typename T>
RK3_Step2_functional3D<T>* RK3_Step2_functional3D<T>::clone() const {
    return new RK3_Step2_functional3D<T>(*this);
}

template<typename T>
void RK3_Step2_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // phi_n
    modified[1] = modif::nothing; // phi_1
    modified[2] = modif::nothing; // phi_1_adv
    modified[3] = modif::staticVariables; // phi_2
}

template<typename T>
BlockDomain::DomainT RK3_Step2_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope; // Everything is local, no communication needed.
}

/* ******** RK3_Step3_functional3D ****************************************** */

template<typename T>
void RK3_Step3_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size() == 4 );
    ScalarField3D<T>& phi_n = *fields[0];
    ScalarField3D<T>& phi_2 = *fields[1];
    ScalarField3D<T>& phi_2_adv = *fields[2];
    ScalarField3D<T>& volfracField_RK = *fields[3];
    Dot3D offset_phi_2     = computeRelativeDisplacement(phi_n, phi_2);
    Dot3D offset_phi_2_adv  = computeRelativeDisplacement(phi_n, phi_2_adv);
    Dot3D offset_volfracField_RK      = computeRelativeDisplacement(phi_n, volfracField_RK);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                volfracField_RK.get(iX+offset_volfracField_RK.x,iY+offset_volfracField_RK.y,iZ+offset_volfracField_RK.z)
                    = 1./3. * phi_n.get(iX,iY,iZ) +
                      2./3. * phi_2.get(iX+offset_phi_2.x,iY+offset_phi_2.y,iZ+offset_phi_2.z) +
                      2./3. * phi_2_adv.get(iX+offset_phi_2_adv.x,iY+offset_phi_2_adv.y,iZ+offset_phi_2_adv.z);
            }
        }
    }
}

template<typename T>
RK3_Step3_functional3D<T>* RK3_Step3_functional3D<T>::clone() const {
    return new RK3_Step3_functional3D<T>(*this);
}

template<typename T>
void RK3_Step3_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // phi_n
    modified[1] = modif::nothing; // phi_2
    modified[2] = modif::nothing; // phi_2_adv
    modified[3] = modif::staticVariables; // volfracField_RK
}

template<typename T>
BlockDomain::DomainT RK3_Step3_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope; // Everything is local, no communication needed.
}



}



#endif  // ADVECTION_DIFFUSION_FD_HH
