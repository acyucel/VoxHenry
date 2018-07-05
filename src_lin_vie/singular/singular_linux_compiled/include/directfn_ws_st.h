/*
 * Copyright (c) 2010 Athanasios Polimeridis
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU LGPL for more details.
 */


#if !defined _DIRECTFN_WS_ST_H_
#define _DIRECTFN_WS_ST_H_

#include <math.h>
#include <complex>
#include "classes.h"

using namespace std;

// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************
//complex<double> directfn_ws_st  ( const double r1[], const double r2[], const double r3[], const int N1, const int N2, const int N3, const int N4, const double k0 );
namespace Triangular_ST
{
	void psi_limits_ws_st ( int argument, double *psi_A, double *psi_B );
	void eta_limits_ws_st ( int argument, double PSI, double *ETA_A, double *ETA_B );

	complex<double> n_functions_ws_st ( int argument, double PSI, double ETA, int kk, Quadrature &quad_LAMBDA, Quadrature &quad_RHO, Geometry_triangle &geom  );

	complex<double> a_functions_ws_st ( double PSI, double ETA, double LAMBDA, int argument, Quadrature &quad_RHO,Geometry_triangle &geom );

	void subtriangles_ws_st(double ETA_sub, double XI_sub, int argument, double *ETA, double *XI);

	void simplex_ws_st(double *xi_1p, double *xi_2p, double *xi_3p, double *xi_1q, double *xi_2q, double *xi_3q, double ETA_p, double XI_p, double ETA_q, double XI_q);
}

complex<double> Kernel(double rp[],double rq[],Geometry_triangle &);




#endif
