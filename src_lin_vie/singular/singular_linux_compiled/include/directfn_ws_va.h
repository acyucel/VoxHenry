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


#if !defined _DIRECTFN_WS_VA_H_
#define _DIRECTFN_WS_VA_H_

#include <math.h>
#include <complex>
#include "classes.h"

using namespace std;


// *********************************************
//			DECLARATION OF FUNCTIONS
// *********************************************


void simplex_ws_va(double *xi_1p, double *xi_2p, double *xi_3p, double *xi_1q, double *xi_2q, double *xi_3q, double Lambda, double Psi, double Theta_p, double Theta_q);

complex<double> Kernel(double rp[], double rq[],Geometry_triangle &);
void gl_quad ( int n, double x[], double w[] );



#endif 
