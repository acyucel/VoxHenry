#include <iostream>
#include "classes.h"
#include "constants.h"
#include "quadric_ws_st.h"
#include <math.h>
using namespace std;

complex<double> quadric_ws_st ( const int N1, const int N2, const int N3, const int N4, Geometry &geom) 
{

        double PSI_A, PSI_B, PSI, J_PSI;
        double U_A, U_B, U, J_U;

        complex<double> I_PSI_const, I_U_const, I_const;

        complex<double> Isub_const[4];

        complex<double> I_subdomain_const[6];

        complex<double> Kernel_const;

        double* w_psi;  // Pointers to arrays of integration points and weights
        double* w_u;
        double* z_psi;
        double* z_u;

        
    // Quadrature parameters

        Quadrature quad_PSI(N1);
		Quadrature quad_U(N2);
		Quadrature quad_LAMBDA(N3);
		Quadrature quad_RHO(N4);

		w_psi = quad_PSI.w;
        z_psi = quad_PSI.z;
        w_u = quad_U.w;
        z_u = quad_U.z;

    for ( int kk = 1 ; kk < 5 ; kk++ )
    {
        for ( int m = 1 ; m < 7 ; m++ )
        {
                     Quadrilateral_ST::psi_limits_ws_st ( m, &PSI_A, &PSI_B );
                     J_PSI = (PSI_B - PSI_A) / double(2.0);
					
                     //
                     I_PSI_const = complex<double>( 0.0 , 0.0 );

                     for ( int n_PSI = 0 ; n_PSI <  N1 ; n_PSI++ )
                     {
                         PSI = ( (PSI_B - PSI_A) / double(2.0) ) * z_psi[n_PSI] + ( (PSI_B + PSI_A) / double(2.0) );
                         //
                         Quadrilateral_ST::u_limits_ws_st ( m, PSI , &U_A, &U_B );
						 
                         J_U = (U_B - U_A) / double(2.0);
                         //

                         I_U_const = complex<double>( 0.0 , 0.0 );

                         for ( int n_U = 0 ; n_U <  N2 ; n_U++ )
                         {
                             U = ( (U_B - U_A) / double(2.0) ) * z_u[n_U] + ( (U_B + U_A) / double(2.0) );
                             //
                             Kernel_const = sin(PSI) * Quadrilateral_ST::n_functions_ws_st ( m, PSI, U, kk, quad_LAMBDA, quad_RHO, geom );
							 
                             //
                             I_U_const = I_U_const + w_u[n_U] * Kernel_const;
                         }
                         I_U_const = J_U * I_U_const;
                         //
                         I_PSI_const = I_PSI_const + w_psi[n_PSI] * I_U_const;
						 
                     }

                     I_subdomain_const[m-1] = J_PSI * I_PSI_const;

                 }

                 for ( int index = 0 ; index < 6 ; index++ )
                 {
                     Isub_const[kk-1] += I_subdomain_const[index];
					 
                 }

             }

         I_const = Isub_const[0] + Isub_const[1] + Isub_const[2] + Isub_const[3];
    
		
         // Final
         return I_const;
}

