
#include <iostream>
#include "classes.h"
#include "QUADRIC_WS_EA.h"
#include <math.h>
using namespace std;

// ***********************************************************************
using namespace Quadrilateral_EA;
// ***********************************************************************
namespace Quadrilateral_EA
{

	complex<double> A_functions_ws_ea (double PSI, double THETA, double U, double LAMBDA_L, Quadrature &quad_LAMBDA,Geometry &geom )
	{

		double* w;    // Pointers to arrays of integration points and weights
		double* z_LAMBDA;

		const int N_LAMBDA = quad_LAMBDA.N;
		w        = quad_LAMBDA.w;
		z_LAMBDA = quad_LAMBDA.z;
    
		//
		complex<double> I_LAMBDA;
		//
		double LAMBDA_a = double(0.0);
		double LAMBDA_b = LAMBDA_L;
		//
		double J_LAMBDA = (LAMBDA_b - LAMBDA_a) / double(2.0);
    //
		double u_p, v_p;
		double u_q, v_q;
		double JpJq;
		double LAMBDA;
		double LAMBDA_Jacobian;

		double rp[3],rq[3];
		complex<double> K;
    //
		I_LAMBDA = complex<double>((double)0.0 ,(double)0.0 );

		for ( int n_LAMBDA = 0 ; n_LAMBDA<  N_LAMBDA ; n_LAMBDA++ )
	    {
			LAMBDA = ( (LAMBDA_b - LAMBDA_a) / double(2.0) ) * z_LAMBDA[n_LAMBDA] + ( (LAMBDA_b + LAMBDA_a) / double(2.0) );
			// Simplex coordinates
			Quadrilateral_EA::Simplex_ws_ea( &u_p, &v_p, &u_q, &v_q, LAMBDA, PSI, THETA, U);
            
            
			JpJq = geom.Jacobians(u_p, v_p, u_q, v_q);
			geom.position_vectors(u_p, v_p, u_q, v_q,rp,rq);

			K = Kernel(rp,rq,geom);
			// lambda Jacobian

			LAMBDA_Jacobian = LAMBDA * LAMBDA;
			//
			I_LAMBDA = I_LAMBDA + w[n_LAMBDA] * LAMBDA_Jacobian * JpJq * K;
		}
		I_LAMBDA = J_LAMBDA * I_LAMBDA;
		// Final output
		return I_LAMBDA;
}
}