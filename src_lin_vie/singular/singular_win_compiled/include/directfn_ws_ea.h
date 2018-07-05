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


#if !defined _DIRECTFN_WS_EA_H_
#define _DIRECTFN_WS_EA_H_

#include <math.h>
#include <complex>
#include "classes.h"

using namespace std;


// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************
namespace Triangular_EA
{
	void psi_limits_ws_ea ( int argument, double theta, double *psi_A, double *psi_B );
	void theta_limits_ws_ea ( int argument, double *theta_A, double *theta_B );

	complex<double> n_functions_ws_ea ( int argument, double PSI, double THETA, int B, Quadrature &, Quadrature &, Geometry_triangle &);

	complex<double> a_functions_ws_ea ( double PSI, double THETA, double ETA, double LAMBDA_L, int B, Quadrature &, Geometry_triangle & );

	void simplex_ws_ea(double *xi_1p, double *xi_2p, double *xi_3p, double *xi_1q, double *xi_2q, double *xi_3q, double Lambda, double Psi, double Theta, double Eta_, int argument );
}

complex<double> Kernel(double rp[],double rq[],Geometry_triangle &);

void gl_quad ( int n, double x[], double w[] );


#endif
