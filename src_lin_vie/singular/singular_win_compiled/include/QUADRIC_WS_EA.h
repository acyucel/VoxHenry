#if !defined _QUADRIC_WS_EA_H_
#define _QUADRIC_WS_EA_H_

#include <math.h>
#include <complex>
#include "classes.h"

using namespace std;


// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************
complex<double> quadric_ws_ea ( const int N1,const int N2, const int N3, const int N4 , Geometry &geom);

complex<double> Kernel(double rp[],double rq[],Geometry &geom);


namespace Quadrilateral_EA
{
	complex<double> N_functions_ws_ea ( int argument, double PSI, double THETA, Quadrature &, Quadrature &, Geometry &);
	complex<double> A_functions_ws_ea ( double PSI, double THETA, double ETA, double LAMBDA_L, Quadrature &, Geometry & );
	void PSI_limits_ws_ea ( int argument, double theta, double *psi_A, double *psi_B );
	void THETA_limits_ws_ea ( int argument, double *theta_A, double *theta_B );
	void Simplex_ws_ea(double *u_p, double *v_p, double *u_q, double *v_q, double Lambda, double Psi, double Theta, double u_ );
}


// ********
//  Macros
// ********

#define round(x)	        ( (int)((x)+ 0.5) )
#define RoundToZero(x)		( fabs(x) > 1.0e-12 ? (x) : 0.0 )

#endif
