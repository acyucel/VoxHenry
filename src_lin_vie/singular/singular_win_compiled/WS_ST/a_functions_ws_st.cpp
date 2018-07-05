#include <iostream>
#include "classes.h"
#include "quadric_ws_st.h"
#include <math.h>
using namespace std;

using namespace Quadrilateral_ST;

namespace Quadrilateral_ST
{
	complex<double> a_functions_ws_st(double PSI, double U, double LAMBDA,int argument, Quadrature &quad_RHO,Geometry &geom)
	{
	double* w;    // Pointers to arrays of integration points and weights
    double* z_RHO;

	const int N_RHO = quad_RHO.N;
	w        = quad_RHO.w;
    z_RHO = quad_RHO.z;

	complex<double> I_RHO;
    //
    double RHO_a = double(0.0);
    double RHO_b = LAMBDA;
    //
    double J_RHO = (RHO_b - RHO_a) / double(2.0);
    //
    double U_p_sub = U;
    double V_p_sub  = LAMBDA * sin(PSI) - double(1.0);
    //
    double U_q_sub, V_q_sub, U_p, V_p, U_q, V_q;
    double RHO;
    double RHO_Jacobian;
	//
	double JpJq;
	double rp[3],rq[3];
    complex<double> K;
    //
    I_RHO = complex<double>((double)0.0 ,(double)0.0 );

    for ( int n_RHO = 0 ; n_RHO<  N_RHO ; n_RHO++ )
    {
        RHO = ( (RHO_b - RHO_a) / double(2.0) ) * z_RHO[n_RHO] + ( (RHO_b + RHO_a) / double(2.0) );
        //
        U_q_sub = RHO  * cos(PSI) + U_p_sub;
        V_q_sub  = -RHO * sin(PSI) + V_p_sub;

        //
        subtriangles_ws_st(U_p_sub, V_p_sub, argument, &U_p, &V_p);
        subtriangles_ws_st(U_q_sub, V_q_sub, argument, &U_q, &V_q);

		// Jacobians
        JpJq = geom.Jacobians(U_p, V_p, U_q, V_q);

		// Position vectors
		geom.position_vectors(U_p, V_p, U_q, V_q,rp,rq);
		
		// Kernel 
		K = Kernel(rp, rq, geom);
		
        // rho Jacobian
        RHO_Jacobian = RHO;
        //
        I_RHO = I_RHO + w[n_RHO] * JpJq * RHO_Jacobian * K;
    }
    I_RHO = J_RHO * I_RHO;
    // Final output
    return I_RHO;
	}
}