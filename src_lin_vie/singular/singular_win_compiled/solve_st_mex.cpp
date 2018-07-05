#include "mex.h"
#include "directfn.h"


/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=17) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "17 inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "1 output required.");   
  
  // INPUT

  //  create a pointer to the input vectors 
  
 
  double* r1 = mxGetPr(prhs[0]);
  double* r2 = mxGetPr(prhs[1]);
  double* r3 = mxGetPr(prhs[2]);
  double* r4 = mxGetPr(prhs[3]);
  
  int N1 = (int)mxGetScalar(prhs[4]);
  int N2 = (int)mxGetScalar(prhs[5]);
  int N3 = (int)mxGetScalar(prhs[6]);
  int N4 = (int)mxGetScalar(prhs[7]);   
  //  
  double k0 = mxGetScalar(prhs[8]);
  double dx = mxGetScalar(prhs[9]);
  double* rq_c = mxGetPr(prhs[10]);
  double* rp_c = mxGetPr(prhs[11]);
  double* nq = mxGetPr(prhs[12]);
  double* np = mxGetPr(prhs[13]);
  
  int kernel_type = (int)mxGetScalar(prhs[14]);
  int integral_type = (int)mxGetScalar(prhs[15]);
  int volume_ker = (int)mxGetScalar(prhs[16]);
 
 
  
  int size;
  size = 1;
  
  
  // OUTPUT
  //  set the output pointer to the output matrix 
  plhs[0] = mxCreateDoubleMatrix(size, 1, mxCOMPLEX);  
  
  //  create a C pointer to a copy of the output matrix 
  double* I_DEr = mxGetPr(plhs[0]);
  double* I_DEi = mxGetPi(plhs[0]);
  
  complex<double> *I_DE = new complex<double>[size];
  /*  call the C subroutine */
  create_ST( r1, r2, r3, r4, N1, N2, N3, N4, k0, dx, rq_c, rp_c, nq, np, kernel_type, integral_type , volume_ker, I_DE);
  
  for (int i = 0; i < size; i++)
  {
      I_DEr[i] = real(I_DE[i]);
      I_DEi[i] = imag(I_DE[i]);
  }
   
  delete I_DE;
}
