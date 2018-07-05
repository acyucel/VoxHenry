#ifndef _CONSTANTS_h
#define _CONSTANTS_h

#include <math.h>
#include <complex>

using namespace std;



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
