

#include <iostream>
#include "classes.h"
#include "constants.h"
#include "QUADRIC_WS_EA.h"

#include <math.h>
using namespace std;

// ***********************************************************************
using namespace Quadrilateral_EA;
// ***********************************************************************
namespace Quadrilateral_EA
{
	complex<double> N_functions_ws_ea (int argument, double PSI, double THETA, Quadrature &quad_U, Quadrature &quad_LAMBDA, Geometry &geom )
	{
		double* w_u;    // Pointers to arrays of integration points and weights
		double* z_u;

		complex<double> I_U, I_U_1, I_U_2, N;
		double U_1, U_2, U;
		double J_U, J_U_1, J_U_2;
		double U_a, U_b , U_1a, U_1b, U_2a, U_2b;

		const int N_u = quad_U.N;
		w_u = quad_U.w;
		z_u = quad_U.z;
     
     double u_1_Psi,u_2_Psi,u_1_theta,u_2_theta;
     double LAMBDA_L1, LAMBDA_L2;
     complex<double> A_L1, A_L2, A_L3;

     switch(argument)
     {
     case 1: //1
         
         //
         I_U = complex<double>((double)0.0 ,(double)0.0 );
         //
         U_a = double(-1.0);
         U_b = double(1.0);

         //
         J_U = ( U_b - U_a) / double(2.0);
        

         for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
         {
             U = (( U_b - U_a) / double(2.0)) * z_u[n_u] + (( U_b + U_a) / double(2.0)) ;
             //
             LAMBDA_L1 = ( double(1.0) + U ) / ( cos (THETA) * cos (PSI) );
             A_L1 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U,LAMBDA_L1,quad_LAMBDA,geom );

             I_U = I_U + w_u[n_u] * A_L1;
         }
         I_U = J_U * I_U;
         //
         N = I_U;
         break;

      case 2: //2,5
		  u_1_Psi = (double(2.0) * cos (THETA)) / tan (PSI) - double(1.0);

          //
          I_U_1 = complex<double>((double)0.0 ,(double)0.0 );
          I_U_2 = complex<double>((double)0.0 ,(double)0.0 );
          //
          U_1a = double(-1.0);
          U_1b = u_1_Psi;

          U_2a = u_1_Psi;
          U_2b = double(1.0);
          //
          J_U_1 = ( U_1b - U_1a) / double(2.0);
          J_U_2 = ( U_2b - U_2a) / double(2.0);

          for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
          {
              U_1 = (( U_1b - U_1a) / double(2.0)) * z_u[n_u] + (( U_1b + U_1a) / double(2.0)) ;
              U_2 = (( U_2b - U_2a) / double(2.0)) * z_u[n_u] + (( U_2b + U_2a) / double(2.0)) ;
              //
              LAMBDA_L1 = ( double(1.0) + U_1 ) / ( cos (THETA) * cos (PSI) );
              A_L1 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_1,LAMBDA_L1,quad_LAMBDA, geom);

              LAMBDA_L2 = ( double(2.0) ) / sin(PSI);
              A_L2 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_2,LAMBDA_L2,quad_LAMBDA, geom);
              //
              I_U_1 = I_U_1 + w_u[n_u] * A_L1;
              I_U_2 = I_U_2 + w_u[n_u] * A_L2;
          }
          I_U_1 = J_U_1 * I_U_1;
          I_U_2 = J_U_2 * I_U_2;
          //
          N = I_U_1 + I_U_2;

          break;

      case 3: // 3,10
          
		  u_1_theta = double(2.0) * (tan ( M_PI/double(2.0) - THETA ) ) - double(1.0);
          //

          I_U_1 = complex<double>((double)0.0 ,(double)0.0 );
          I_U_2 = complex<double>((double)0.0 ,(double)0.0 );
          //
          U_1a = double(-1.0);
          U_1b = u_1_theta;

          U_2a = u_1_theta;
          U_2b = double(1.0);
          //
          J_U_1 = ( U_1b - U_1a) / double(2.0);
          J_U_2 = ( U_2b - U_2a) / double(2.0);

          for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
          {
              U_1 = (( U_1b - U_1a) / double(2.0)) * z_u[n_u] + (( U_1b + U_1a) / double(2.0)) ;
              U_2 = (( U_2b - U_2a) / double(2.0)) * z_u[n_u] + (( U_2b + U_2a) / double(2.0)) ;
              //
              LAMBDA_L1 = ( double(1.0) + U_1 ) / ( cos (THETA) * cos (PSI) );
              A_L1 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_1,LAMBDA_L1,quad_LAMBDA, geom);

              LAMBDA_L2 = ( double(2.0) ) / (sin(THETA) * cos (PSI));
              A_L2 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_2,LAMBDA_L2,quad_LAMBDA, geom);
              //
              I_U_1 = I_U_1 + w_u[n_u] * A_L1;
              I_U_2 = I_U_2 + w_u[n_u] * A_L2;
          }
          I_U_1 = J_U_1 * I_U_1;
          I_U_2 = J_U_2 * I_U_2;
          //
          N = I_U_1 + I_U_2;

          break;

      case 4: // 4,6+9

          u_1_theta = double(2.0) * (tan ( M_PI/double(2.0) - THETA ) ) - double(1.0);
		  u_1_Psi = (double(2.0) * cos (THETA)) / tan (PSI) - double(1.0);
          //

          I_U_1 = complex<double>((double)0.0 ,(double)0.0 );
          I_U_2 = complex<double>((double)0.0 ,(double)0.0 );

          //
          U_1a = double(-1.0);
          U_1b = u_1_Psi;

          U_2a = u_1_Psi;
          U_2b = double(1.0);

          //
          J_U_1 = ( U_1b - U_1a) / double(2.0);
          J_U_2 = ( U_2b - U_2a) / double(2.0);

          for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
          {
              U_1 = (( U_1b - U_1a) / double(2.0)) * z_u[n_u] + (( U_1b + U_1a) / double(2.0)) ;
              U_2 = (( U_2b - U_2a) / double(2.0)) * z_u[n_u] + (( U_2b + U_2a) / double(2.0)) ;
              //
              LAMBDA_L1 = ( double(1.0) + U_1 ) / ( cos (THETA) * cos (PSI) );
              A_L1 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_1,LAMBDA_L1,quad_LAMBDA, geom);

              LAMBDA_L2 = ( double(2.0) ) / (sin(PSI));
              A_L2 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_2,LAMBDA_L2,quad_LAMBDA, geom);
              //
              I_U_1 = I_U_1 + w_u[n_u] * A_L1;
              I_U_2 = I_U_2 + w_u[n_u] * A_L2;
          }
          I_U_1 = J_U_1 * I_U_1;
          I_U_2 = J_U_2 * I_U_2;
          //
          N = I_U_1 + I_U_2;

          break;

      case 5: // 8,11
          u_2_theta = double(2.0) * (tan ( M_PI/double(2.0) - THETA ) ) + double(1.0);
          //

          I_U_1 = complex<double>((double)0.0 ,(double)0.0 );
          I_U_2 = complex<double>((double)0.0 ,(double)0.0 );
          //
          U_1a = -double(1.0);
          U_1b = u_2_theta;

          U_2a = u_2_theta;
          U_2b = double(1.0);
          //
          J_U_1 = ( U_1b - U_1a) / double(2.0);
          J_U_2 = ( U_2b - U_2a) / double(2.0);

          for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
          {
              U_1 = (( U_1b - U_1a) / double(2.0)) * z_u[n_u] + (( U_1b + U_1a) / double(2.0)) ;
              U_2 = (( U_2b - U_2a) / double(2.0)) * z_u[n_u] + (( U_2b + U_2a) / double(2.0)) ;
              //
              LAMBDA_L1 = ( double(2.0) ) / (sin(THETA) * cos(PSI));
              A_L1 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_1,LAMBDA_L1,quad_LAMBDA, geom);

              LAMBDA_L2 = ( U_2 - double(1.0)) / (cos(THETA) * cos(PSI));
              A_L2 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_2,LAMBDA_L2,quad_LAMBDA, geom);
              //
              I_U_1 = I_U_1 + w_u[n_u] * A_L1;
              I_U_2 = I_U_2 + w_u[n_u] * A_L2;
          }
          I_U_1 = J_U_1 * I_U_1;
          I_U_2 = J_U_2 * I_U_2;
          //
          N = I_U_1 + I_U_2;

          break;


	  case 6:  // 10+15, 12
 
		  u_2_Psi = (double(2.0) * cos (THETA)) / tan (PSI)  + double(1.0);
          //

          I_U_1 = complex<double>((double)0.0 ,(double)0.0 );
          I_U_2 = complex<double>((double)0.0 ,(double)0.0 );
          //
          U_1a = double(-1.0);
          U_1b = u_2_Psi;

          U_2a = u_2_Psi;
          U_2b = double(1.0);

          //
          J_U_1 = ( U_1b - U_1a) / double(2.0);
          J_U_2 = ( U_2b - U_2a) / double(2.0);
 

          for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
          {
              U_1 = (( U_1b - U_1a) / double(2.0)) * z_u[n_u] + (( U_1b + U_1a) / double(2.0)) ;
              U_2 = (( U_2b - U_2a) / double(2.0)) * z_u[n_u] + (( U_2b + U_2a) / double(2.0)) ;
              //
              LAMBDA_L1 = ( double(2.0) ) / (sin(PSI));
              A_L1 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_1,LAMBDA_L1,quad_LAMBDA, geom);

              LAMBDA_L2 = ( U_2 - double(1.0)) / (cos(THETA) * cos(PSI));
              A_L2 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U_2,LAMBDA_L2,quad_LAMBDA, geom);
              //

              I_U_1 = I_U_1 + w_u[n_u] * A_L1;
              I_U_2 = I_U_2 + w_u[n_u] * A_L2;
          }
          I_U_1 = J_U_1 * I_U_1;
          I_U_2 = J_U_2 * I_U_2;
          //
          N = I_U_1 + I_U_2;

          break;

	  case 7: // 13
		 //


         I_U = complex<double>((double)0.0 ,(double)0.0 );
         //
		 U_a = double(-1.0);
         U_b = double(1.0);
         

         //
         J_U = ( U_b - U_a) / double(2.0);
        

         for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
         {
             U = (( U_b - U_a) / double(2.0)) * z_u[n_u] + (( U_b + U_a) / double(2.0)) ;
             //
             LAMBDA_L2 = ( U - double(1.0)) / (cos(THETA) * cos(PSI));
             A_L2 = Quadrilateral_EA::A_functions_ws_ea(PSI,THETA,U,LAMBDA_L2,quad_LAMBDA, geom);

             I_U = I_U + w_u[n_u] * A_L2;
         }
         I_U = J_U * I_U;
         //
         N = I_U;
         break;

	  case 8:  // 14,16

          u_2_Psi = (double(2.0) * cos (THETA)) / tan (PSI)  + double(1.0);
          //

          I_U_1 = complex<double>((double)0.0 ,(double)0.0 );
          I_U_2 = complex<double>((double)0.0 ,(double)0.0 );
          //
          U_1a = u_2_Psi;
          U_1b = double(1.0);

          U_2a = double(-1.0);
          U_2b = u_2_Psi;
          //
          J_U_1 = ( U_1b - U_1a) / double(2.0);

          J_U_2 = ( U_2b - U_2a) / double(2.0);

          for ( int n_u = 0 ; n_u <  N_u ; n_u++ )
          {
              U_1 = (( U_1b - U_1a) / double(2.0)) * z_u[n_u] + (( U_1b + U_1a) / double(2.0)) ;
              U_2 = (( U_2b - U_2a) / double(2.0)) * z_u[n_u] + (( U_2b + U_2a) / double(2.0)) ;
              //
              LAMBDA_L1 = ( U_1 - double(1.0)) / (cos(THETA) * cos(PSI));
              A_L1 = Quadrilateral_EA::A_functions_ws_ea (PSI,THETA,U_1,LAMBDA_L1,quad_LAMBDA, geom);

              LAMBDA_L2 = ( double(2.0) ) / (sin(PSI));
              A_L2 = Quadrilateral_EA::A_functions_ws_ea (PSI,THETA,U_2,LAMBDA_L2,quad_LAMBDA, geom);
              //
              I_U_1 = I_U_1 + w_u[n_u] * A_L1;
              I_U_2 = I_U_2 + w_u[n_u] * A_L2;
          }
          I_U_1 = J_U_1 * I_U_1;
          I_U_2 = J_U_2 * I_U_2;
          //
          N = I_U_1 + I_U_2;

          break;

     }

     // Final output
     
     return N;
	 }
}

