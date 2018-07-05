#include <iostream>
#include "constants.h"
#include "classes.h"
#include "quadric_ws_st.h"
#include <math.h>
using namespace std;

using namespace Quadrilateral_ST;

namespace Quadrilateral_ST
{
	complex<double> n_functions_ws_st  (int argument, double PSI, double U, int kk, Quadrature &quad_LAMBDA, Quadrature &quad_RHO, Geometry &geom )
	{
    double* w;    // Pointers to arrays of integration points and weights
    double* z_LAMBDA;

    complex<double> A_L1;
    complex<double> J_LAMBDA, I_LAMBDA;
    double LAMBDA_a, LAMBDA_b, LAMBDA;

    const int N_LAMBDA = quad_LAMBDA.N;
    w = quad_LAMBDA.w;
    z_LAMBDA = quad_LAMBDA.z;

     
     I_LAMBDA = complex<double>((double)0.0 ,(double)0.0 );
         //
     LAMBDA_a = double(0.0);
     lambda_limits_ws_st (argument,PSI, U, &LAMBDA_b );

         //
     J_LAMBDA = ( LAMBDA_b - LAMBDA_a) / double(2.0);
     for ( int n_LAMBDA = 0 ; n_LAMBDA <  N_LAMBDA ; n_LAMBDA++ )
     {
         LAMBDA = (( LAMBDA_b - LAMBDA_a) / double(2.0)) * z_LAMBDA[n_LAMBDA] + (( LAMBDA_b + LAMBDA_a) / double(2.0)) ;
          
         //
         A_L1 = a_functions_ws_st(PSI,U,LAMBDA,kk,quad_RHO, geom);
         //
         I_LAMBDA = I_LAMBDA + w[n_LAMBDA] * A_L1;
     }
     I_LAMBDA = J_LAMBDA * I_LAMBDA;
     // Final output
     return I_LAMBDA;
	}
}

