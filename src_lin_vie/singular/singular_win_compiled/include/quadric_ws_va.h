#if !defined _QUADRIC_WS_VA_H_
#define _QUADRIC_WS_VA_H_

#include <complex>
#include "classes.h"

using namespace std;


// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************
complex<double> quadric_ws_va ( const int N1, const int N2, const int N3, const int N4,Geometry &geom);

namespace Quadrilateral_VA
{
	void theta_p_limits ( int argument, double *theta_A, double *theta_B );
	void theta_q_limits ( int argument, double *theta_A, double *theta_B );
	void rho_p_limit(int argument,double theta_p, double *L_p);
	void rho_q_limit(int argument,double theta_q, double *L_q);
}
complex<double> Kernel(double rp[],double rq[],Geometry &geom);
//complex<double> Kernel_linear_voxel  (int lp, int lq, double rp[], double rq[], Geometry &geom);
#endif