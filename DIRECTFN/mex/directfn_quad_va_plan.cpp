#include "mex.h"
#include "directfn_quad.h"
#include <memory>
#include <iostream>
#include <iomanip>
#include <string>

using  std::string;
using namespace Directfn;

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=13) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "13 inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "1 output required.");   
  
  // INPUT

  //  create a pointer to the input vectors 
  
  double* r1 = mxGetPr(prhs[0]);
  double* r2 = mxGetPr(prhs[1]);
  double* r3 = mxGetPr(prhs[2]);
  double* r4 = mxGetPr(prhs[3]);
  double* r5 = mxGetPr(prhs[4]);
  double* r6 = mxGetPr(prhs[5]);
  double* r7 = mxGetPr(prhs[6]);

  int N1 = (int)mxGetScalar(prhs[7]);  
  int N2 = (int)mxGetScalar(prhs[8]);
  int N3 = (int)mxGetScalar(prhs[9]);
  int N4 = (int)mxGetScalar(prhs[10]);

  double k0 = mxGetScalar(prhs[11]);

  const string btt_type = mxArrayToString(prhs[12]);
  int len;
  int status = -1;
  unique_ptr<dcomplex[]> I_DE(nullptr);

  if (btt_type == string("Constant")) {
	  len = QuadrilateralKernel_PlanarScalar::constexpr_size();
	  I_DE.reset(new dcomplex[len]);
	  status = directfn_quad_va_plan<QuadrilateralKernel_PlanarScalar>(r1, r2, r3, r4, r5, r6, r7, N1, N2, N3, N4, k0, I_DE.get());
  }
  else if (btt_type == string("Vector_WS")) {
	  len = QuadrilateralKernel_PlanarVectorWS::constexpr_size();
	  I_DE.reset(new dcomplex[len]);
	  status = directfn_quad_va_plan<QuadrilateralKernel_PlanarVectorWS>(r1, r2, r3, r4, r5, r6, r7, N1, N2, N3, N4, k0, I_DE.get());
  }
  else if (btt_type == string("Vector_SS")) {
	  len = QuadrilateralKernel_PlanarVectorSS::constexpr_size();
	  I_DE.reset(new dcomplex[len]);
	  status = directfn_quad_va_plan<QuadrilateralKernel_PlanarVectorSS>(r1, r2, r3, r4, r5, r6, r7, N1, N2, N3, N4, k0, I_DE.get());
  }
  else mexErrMsgTxt("Wrong basis/testing parameter!");
  
  
//   mexPrintf("z_psi =  (%f,%f,%f) \n", z_psi[7],z_psi[8],z_psi[9]);
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(len, 1, mxCOMPLEX);  
  
  //  create a C pointer to a copy of the output matrix 
  double* I_DEr = mxGetPr(plhs[0]);
  double* I_DEi = mxGetPi(plhs[0]);
  
  //mexPrintf("%d %d",real(I_DE),imag(I_DE));
   for (int i = 0; i < len; i++)
  {
      I_DEr[i] = real(I_DE[i]);
      I_DEi[i] = imag(I_DE[i]);
  }
  
}
