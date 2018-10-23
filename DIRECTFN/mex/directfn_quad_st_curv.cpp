#include "mex.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "directfn_quad.h"

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
  if(nrhs!=15) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "15 inputs required.");
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
  double* r8 = mxGetPr(prhs[7]);
  double* r9 = mxGetPr(prhs[8]);
 
    
  int N1 = (int)mxGetScalar(prhs[9]);
  int N2 = (int)mxGetScalar(prhs[10]);
  int N3 = (int)mxGetScalar(prhs[11]);
  int N4 = (int)mxGetScalar(prhs[12]);
  //  get the scalar input ko 
  double k0 = mxGetScalar(prhs[13]);

  
  const string btt_type = mxArrayToString(prhs[14]);
  int len;
  int status = -1;
  unique_ptr<dcomplex[]> I_DE(nullptr);
  
  if (btt_type == string("Constant")) {
	  len = QuadrilateralKernel_CurvilinearScalar::constexpr_size();
	  I_DE.reset(new dcomplex[len]);
	  status = directfn_quad_st_curv<QuadrilateralKernel_CurvilinearScalar>(r1, r2, r3, r4, r5, r6, r7, r8, r9, N1, N2, N3, N4, k0, I_DE.get());
  }
  else if (btt_type == string("Vector_WS")) {
	  len = QuadrilateralKernel_CurvilinearVectorWS::constexpr_size();
	  I_DE.reset(new dcomplex[len]);
	  status = directfn_quad_st_curv<QuadrilateralKernel_CurvilinearVectorWS>(r1, r2, r3, r4, r5, r6, r7, r8, r9, N1, N2, N3, N4, k0, I_DE.get());
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
