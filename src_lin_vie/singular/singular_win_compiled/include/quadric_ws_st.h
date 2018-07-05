#if !defined _QUADRIC_WS_ST_H_
#define _QUADRIC_WS_EA_H_

#include <complex>
#include "classes.h"

using namespace std;

// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************


complex<double> quadric_ws_st ( const int N1,const int N2,const int N3,const int N4,Geometry &geom);

namespace Quadrilateral_ST
{
	void u_limits_ws_st ( int argument, double PSI , double *U_A, double *U_B );
	void psi_limits_ws_st ( int argument, double *psi_A, double *psi_B );
	void lambda_limits_ws_st ( int argument, double PSI, double U, double *LAMBDA_L );
	complex<double> a_functions_ws_st(double PSI, double U, double LAMBDA,int argument, Quadrature &quad_RHO,Geometry &geom);
	complex<double> n_functions_ws_st  (int argument, double PSI, double U, int kk, Quadrature &quad_LAMBDA, Quadrature &quad_RHO, Geometry &geom );
	void subtriangles_ws_st(double U_sub, double V_sub, int argument, double *U, double *V);
}

complex<double> Kernel(double rp[],double rq[],Geometry &geom);
//complex<double> Kernel_linear_voxel  (int lp, int lq, double rp[], double rq[], Geometry &geom);

#endif