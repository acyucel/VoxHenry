
//

#ifndef _DIRECTFN_CONST_h
#define _DIRECTFN_CONST_h

#include "classes.h"

using namespace std;

// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************

complex<double> quadric_ws_st  ( const int N1, const int N2, const int N3, const int N4, Geometry &geom);

complex<double> quadric_ws_ea ( const int N1,const int N2, const int N3, const int N4, Geometry &geom);

complex<double>  quadric_ws_va  (const int N1,const int N2,const int N3,const int N4, Geometry &geom);

complex<double> directfn_ws_st  ( const double r1[], const double r2[], const double r3[], const int N1, const int N2, const int N3, const int N4, const double k0 );

complex<double> directfn_ws_ea ( const double r1[], const double r2[], const double r3[], const double r4[], const int N1, const int N2, const int N3, const int N4, const double k0 );

complex<double>  directfn_ws_va ( const double r1[], const double r2[], const double r3[], const double r4[],const double r5[],const int N1, const int N2, const int N3, const int N4, const double k0);
void create_ST(const double r1[], const double r2[], const double r3[], const double r4[], int N1 , int N2, int N3,  int N4,double k0, double dx, double rq_c[],double rp_c[],double nq[], double np[], int ker_type,int integral_type, int volume_ker, complex<double> I[]);

void create_EA(const double r1[],const double r2[],const double r3[],const double r4[],const double r5[],const double r6[],const int N1,const int N2, const int N3, const int N4, double k0, double dx, double rq_c[],double rp_c[],double nq[], double np[], int ker_type, int integral_type, int volume_ker, complex<double> I[]);


void create_VA(const double r1[],const double r2[],const double r3[],const double r4[],const double r5[],const double r6[], const double r7[],const int N1,const int N2, const int N3, const int N4, double k0, double dx, double rq_c[],double rp_c[],double nq[], double np[], int ker_type, int integral_type, int volume_ker, complex<double> I[]);
// **************************************
//			Inline functions
// **************************************

/*
inline
double vector_dot(double x[], double y[])
{
    return x[0]*y[0]+x[1]*y[1]+x[2]*y[2];
}
//
inline
void vector_cross(double x[], double y[], double z[])
{
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
}
*/
// **************************************
//			Mathematical Constants
// **************************************

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830961566084582     /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

// **************************************
//			Physical Constants
// **************************************

const double eo	     =   8.85400e-12;			// free space electric permitivity
const double mo	     =   4.0 * M_PI * 1.0e-7;		// free space magnetic permeability
const double co             =   299792458;			// free spave light velocity
const double Zo             =   376.734;			// free space inpedance



const std::complex<double> Iunit = std::complex<double>( 0.0 , 1.0 );


const double freq = double(75.0)*1.0e+6;

const double v_light_0 = co;
const double mu_0 = mo;
const double epsilon_0 = double(1.0) / (v_light_0 * v_light_0 * mu_0);

const double omega = double(2.0) * M_PI * freq;
//const double ko   = omega /v_light_0;
//const double ko = double(2.0)*M_PI;
#endif
