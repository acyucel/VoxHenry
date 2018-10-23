#include "mex.h"
#include "directfn_triag.h"
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
  if(nrhs!=9) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "9 inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "1 output required.");   
  
  // INPUT

  //  create a pointer to the input vectors 
  
  double* r1 = mxGetPr(prhs[0]);
  double* r2 = mxGetPr(prhs[1]);
  double* r3 = mxGetPr(prhs[2]);
 
    
  int N1 = (int)mxGetScalar(prhs[3]);
  int N2 = (int)mxGetScalar(prhs[4]);
  int N3 = (int)mxGetScalar(prhs[5]);
  int N4 = (int)mxGetScalar(prhs[6]);

  double k0 = mxGetScalar(prhs[7]);

  const string btt_type = mxArrayToString(prhs[8]);
  int status = -1;
  int len;
  unique_ptr<dcomplex[]> I_DE(nullptr);
  
  if (btt_type == string("Constant")) {
	  len = TriangularKernel_Constant_ST::constexpr_size();
	  I_DE.reset(new dcomplex[len]);
	  status = directfn_tri_st_plan<TriangularKernel_Constant_ST>(r1, r2, r3, N1, N2, N3, N4, k0, I_DE.get());
  }
  else if (btt_type == string("RWG_WS")) {
	  len = TriangularKernel_RWG_WS::constexpr_size();
	  I_DE.reset(new dcomplex[len]);
	  status = directfn_tri_st_plan<TriangularKernel_RWG_WS>(r1, r2, r3, N1, N2, N3, N4, k0, I_DE.get());
  }
  else mexErrMsgTxt("Wrong basis/testing parameter!");
    
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(len, 1, mxCOMPLEX);  
  
  //  create a C pointer to a copy of the output matrix 
  double* I_DEr = mxGetPr(plhs[0]);
  double* I_DEi = mxGetPi(plhs[0]);
 

  for (int i = 0; i < len; i++)
  {
      I_DEr[i] = real(I_DE[i]);
      I_DEi[i] = imag(I_DE[i]);
  }
  
}
