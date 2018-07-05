
#include <iostream>
//#include "directfn_inline.h"
#include "classes.h"
#include "constants.h"
#include "QUADRIC_WS_EA.h"
#include <math.h>
using namespace std;

using namespace Quadrilateral_EA;
// ***********************************************************************

// ***********************************************************************
	complex<double> quadric_ws_ea (const int N1,const int N2, const int N3, const int N4,Geometry &geom)
	{

    // ************************************************
	//			DECLARATION OF KEY VARIABLES
	// ************************************************

	// 1. Various

	    complex<double> I_DE;

        double psi_A, psi_B, PSI, J_psi;
        double theta_A, theta_B, THETA;

        double J_theta = double(0.0);

        complex<double> I_psi_const, I_theta_const, I_const;

		complex<double> X;

        complex<double> I_subdomain_const[14];

        complex<double> Kernel_const;

		
	// 3. Quadrature parameters

        Quadrature quad_THETA(N1);
		Quadrature quad_PSI(N2);
		Quadrature quad_U(N3);
		Quadrature quad_LAMBDA(N4);

		
        double* w_psi;  // Pointers to arrays of integration points and weights
        double* w_theta;
        double* z_psi;
        double* z_theta;

         w_psi = quad_PSI.w;
         z_psi = quad_PSI.z;
         w_theta = quad_THETA.w;
         z_theta = quad_THETA.z;
	
	 // for loops
         
         for ( int m = 1 ; m < 9 ; m++ )
            {
             I_theta_const   = complex<double>( 0.0 , 0.0 );

             for ( int n_theta = 0 ; n_theta <  N1 ; n_theta++ )
                {
                 THETA_limits_ws_ea ( m, &theta_A, &theta_B );
                 THETA = ( (theta_B - theta_A) / double(2.0) ) * z_theta[n_theta] + ( (theta_B + theta_A) / double(2.0) );
                 J_theta = (theta_B - theta_A) / double(2.0);
                 //
                 I_psi_const   = complex<double>( 0.0 , 0.0 );

                 for ( int n_psi = 0 ; n_psi <  N2 ; n_psi++ )
                    {
                     PSI_limits_ws_ea ( m, THETA, &psi_A, &psi_B );
                     PSI = ( (psi_B - psi_A) / double(2.0) ) * z_psi[n_psi] + ( (psi_B + psi_A) / double(2.0) );
                     J_psi = (psi_B - psi_A) / double(2.0);

					 X = N_functions_ws_ea (m, PSI, THETA, quad_U, quad_LAMBDA, geom );
					
					 Kernel_const   = cos(PSI) * X;
					
                     I_psi_const   = I_psi_const   + w_psi[n_psi] * Kernel_const;
					 
                 }
                 I_psi_const   = J_psi * I_psi_const;
				 
                 I_theta_const   = I_theta_const   + w_theta[n_theta] * I_psi_const;

             } 

             I_subdomain_const[m-1]   = J_theta * I_theta_const;
             

         }

         for ( int index = 0 ; index < 8 ; index++ )
            {
              I_DE += I_subdomain_const[index];
            }
         // FINAL OUTPUT
		 
		 return I_DE;

}