
#include <iostream>
//#include "directfn_inline.h"
#include "classes.h"
#include <complex>
#include <math.h>
#include "directfn.h"
using namespace std;

void create_ST(const double r1[], const double r2[], const double r3[], const double r4[], int N1, int N2,int N3,int N4,double k0, double dx, double rq_c[],double rp_c[],double nq[], double np[], int ker_type, int integral_type, int volume_ker, complex<double> I[])
{
	Geometry geom;
	geom.ST(r1,r2,r3,r4);
	geom.set_wavenumber(k0);
	geom.set_delta(dx);
	geom.set_centers(rq_c,rp_c);
	geom.set_normales(nq,np);
	geom.set_kerneltype(ker_type);
	geom.set_lp(integral_type);
	geom.set_volume_kernel(volume_ker);

	I[0] = quadric_ws_st(N1, N2, N3, N4, geom);

	
}


