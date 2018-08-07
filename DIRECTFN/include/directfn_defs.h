#ifndef _DIRECTFN_DEFINITIONS_H_
#define _DIRECTFN_DEFINITIONS_H_

#include <complex>
//#include <cfloat>

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

// Defining constants
#ifndef M_PI
const double M_PI = 3.14159265358979323846264338328;
#endif

using dcomplex = std::complex<double>;
const dcomplex Iunit(0.0, 1.0);


///////////////////////////////////////////////////////////////////////////////

}  // End of the

#endif   // _DIRECTFN_DEFINITIONS_H_

// End of the file



//enum class Singularity : int {ST, EA, VA, NS};

///*! RWG_WS - RWG functions, Weakly Singular kernels, 1/R,
// *  RWG_SS - RWG functions, Strongly Singular kernels, 1/R^2   */
//enum class BasisTestingType : int {Constant,
//                                   RWG_WS, RWG_SS, nxRWG_SS,
//                                   RoofTop_WS, RoofTop_SS,
//                                   Linear_Kop, Linear_Nop};

//                                   //Lin_0, Lin_1, Lin_2, Lin_3, Lin_4,
//                                   //Lin_5, Lin_6, Lin_7, Lin_8};

