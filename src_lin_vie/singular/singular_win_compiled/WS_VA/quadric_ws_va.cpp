#include <iostream>
#include "classes.h"
#include "constants.h"
#include "quadric_ws_va.h"
#include <math.h>
using namespace std;

using namespace Quadrilateral_VA;

complex<double> quadric_ws_va ( const int N1,const int N2,const int N3,const int N4, Geometry &geom)
{

        double theta_p_A, theta_p_B, THETA_p, J_theta_p, Lp;
        double theta_q_A, theta_q_B, THETA_q, J_theta_q, Lq;
        double psi_1_A, psi_1_B, PSI_1, J_psi_1;
        double psi_2_A, psi_2_B, PSI_2, J_psi_2;
        double Lambda_1_A, Lambda_1_B, Lambda_1, J_Lambda_1;
        double Lambda_2_A, Lambda_2_B, Lambda_2, J_Lambda_2;
		double U_p,V_p,U_q,V_q;
		double JpJq;
		double rp[3], rq[3];


        complex<double> Jacobian_1, Jacobian_2;
        complex<double> Kernel1_const, F1_const, Kernel2_const, F2_const;

        complex<double> I_theta_p_const, I_theta_q_const, I_psi_const, I_Lambda1_const, I_Lambda2_const;
		complex<double> I_subdomain_const[4];
        complex<double> I;

	
        double* w_theta_p;
        double* z_theta_p;
        double* w_theta_q;
        double* z_theta_q;
        double* w_psi;  // Pointers to arrays of integration points and weights
        double* z_psi;
        double* w_Lambda;
        double* z_Lambda;
	
		
      
       //  Quadrature parameters 
	   /*
		Quadrature quad_THETA_P(N1);
		Quadrature quad_THETA_Q(N2);
		Quadrature quad_PSI(N3);
		Quadrature quad_LAMBDA(N4);
        
        w_theta_p = quad_THETA_P.w;
        z_theta_p = quad_THETA_P.z;
         //
        w_theta_q = quad_THETA_Q.w;
        z_theta_q = quad_THETA_Q.z;
         //
         w_psi = quad_PSI.w;
         z_psi = quad_PSI.z;
         //
         w_Lambda = quad_LAMBDA.w;
         z_Lambda = quad_LAMBDA.z;
	*/
 // 1. Allocate space for the arrays
         w_theta_p = new double [N1];
         z_theta_p = new double [N1];
         //
         w_theta_q = new double [N2];
         z_theta_q = new double [N2];
         //
         w_psi = new double [N3];
         z_psi = new double [N3];
         //
         w_Lambda = new double [N4];
         z_Lambda = new double [N4];
	// 2. Get the weights and abscissas
         GL_1D ( N1, z_theta_p, w_theta_p );
         GL_1D ( N2, z_theta_q, w_theta_q );
         GL_1D ( N3, z_psi, w_psi );
         GL_1D ( N4, z_Lambda, w_Lambda );


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // for loops

		 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////



for ( int m = 1; m < 5; m++)
{
    I_theta_p_const = 0.0 ;
    for ( int n_theta_p = 0 ; n_theta_p < N1 ; n_theta_p++ )
	 {
             Quadrilateral_VA::theta_p_limits(m, &theta_p_A, &theta_p_B);
             //
             THETA_p = ( (theta_p_B - theta_p_A) / double(2.0) ) * z_theta_p[n_theta_p] + ( (theta_p_B + theta_p_A) / double(2.0) );
             J_theta_p = (theta_p_B - theta_p_A) / double(2.0);
             //
             Quadrilateral_VA::rho_p_limit(m, THETA_p, &Lp);
             //

             I_theta_q_const = 0.0;

             for ( int n_theta_q = 0 ; n_theta_q < N2 ; n_theta_q++ )
             {
                 Quadrilateral_VA::theta_q_limits(m, &theta_q_A, &theta_q_B);
                 //
                 THETA_q = ( (theta_q_B - theta_q_A) / double(2.0) ) * z_theta_q[n_theta_q] + ( (theta_q_B + theta_q_A) / double(2.0) );
                 J_theta_q = (theta_q_B - theta_q_A) / double(2.0);
                 //
                 Quadrilateral_VA::rho_q_limit(m, THETA_q, &Lq);
                 //

                 I_psi_const = 0.0;


                 for ( int n_psi = 0 ; n_psi < N3 ; n_psi++ )
                 {
                     //   0 =< PSI <= PSI_1
                     psi_1_A = double(0.0);
                     psi_1_B = atan(Lq/Lp);
                     //
                     PSI_1   = ( (psi_1_B - psi_1_A) / double(2.0) ) * z_psi[n_psi] + ( (psi_1_B + psi_1_A) / double(2.0) );
                     J_psi_1 = (psi_1_B - psi_1_A) / double(2.0);
                     //

                     I_Lambda1_const = 0.0;


                     for ( int n_Lambda1 = 0 ; n_Lambda1 < N4 ; n_Lambda1++ )
                     {
                         Lambda_1_A = double(0.0);
                         Lambda_1_B = Lp / cos(PSI_1);
                         //
                         Lambda_1   = ( (Lambda_1_B - Lambda_1_A) / double(2.0) ) * z_Lambda[n_Lambda1] + ( (Lambda_1_B + Lambda_1_A) / double(2.0) );
                         J_Lambda_1 = (Lambda_1_B - Lambda_1_A) / double(2.0);

                         U_p = Lambda_1*cos(PSI_1)*cos(THETA_p) - double(1.0);
						 V_p = Lambda_1*cos(PSI_1)*sin(THETA_p) - double(1.0);

						 U_q = Lambda_1*sin(PSI_1)*cos(THETA_q) - double(1.0);
						 V_q = Lambda_1*sin(PSI_1)*sin(THETA_q) - double(1.0);

						 // Jacobians of the parametric transformations
						 JpJq = geom.Jacobians(U_p, V_p, U_q, V_q);

						 // Position vectors
		                 geom.position_vectors(U_p, V_p, U_q, V_q,rp,rq);
						 
						 // Kernel
                         Kernel1_const = Kernel(rp,rq,geom);


                         Jacobian_1 = pow(Lambda_1,3.0) * sin(PSI_1) * cos(PSI_1);

                         F1_const = Kernel1_const * Jacobian_1;
                         //
                         I_Lambda1_const = I_Lambda1_const + w_Lambda[n_Lambda1] * JpJq * F1_const;
                     } 
                     I_Lambda1_const = J_Lambda_1 * I_Lambda1_const;


                     //   PSI_1 =< PSI <= pi/2
                     psi_2_A = atan(Lq/Lp);
                     psi_2_B = M_PI / double(2.0);
                     //
                     PSI_2   = ( (psi_2_B - psi_2_A) / double(2.0) ) * z_psi[n_psi] + ( (psi_2_B + psi_2_A) / double(2.0) );
                     J_psi_2 = (psi_2_B - psi_2_A) / double(2.0);
                     //

                     I_Lambda2_const = 0.0;

                     for ( int n_Lambda2 = 0 ; n_Lambda2 < N4 ; n_Lambda2++ )
                     {
                         Lambda_2_A = double(0.0);
                         Lambda_2_B = Lq / sin(PSI_2);
                         //
                         Lambda_2   = ( (Lambda_2_B - Lambda_2_A) / double(2.0) ) * z_Lambda[n_Lambda2] + ( (Lambda_2_B + Lambda_2_A) / double(2.0) );
                         J_Lambda_2 = (Lambda_2_B - Lambda_2_A) / double(2.0);

                         U_p = Lambda_2*cos(PSI_2)*cos(THETA_p) - double(1.0);
						 V_p = Lambda_2*cos(PSI_2)*sin(THETA_p) - double(1.0);

						 U_q = Lambda_2*sin(PSI_2)*cos(THETA_q) - double(1.0);
						 V_q = Lambda_2*sin(PSI_2)*sin(THETA_q) - double(1.0);

						 // Jacobians of the parametric transformations
						 JpJq = geom.Jacobians(U_p, V_p, U_q, V_q);

						 // Distance function
		                 geom.position_vectors(U_p, V_p, U_q, V_q,rp,rq);
						 
						 // Kernel
                         Kernel2_const = Kernel(rp,rq,geom);

                         Jacobian_2 = pow(Lambda_2,3.0) * sin(PSI_2) * cos(PSI_2);

                         F2_const = Kernel2_const * Jacobian_2;

                         //
                         I_Lambda2_const = I_Lambda2_const + w_Lambda[n_Lambda2] * JpJq * F2_const;

                     } // for ( int n_Lambda1 = 0 ; Lambda1 < N4 ; Lambda1++ )
                     I_Lambda2_const = J_Lambda_2 * I_Lambda2_const;
                     //
                     I_psi_const = I_psi_const + w_psi[n_psi] * ( J_psi_1 * I_Lambda1_const + J_psi_2 * I_Lambda2_const);
                 } // for ( int n_psi = 0 ; n_psi < N3; n_psi++ )
                 I_theta_q_const = I_theta_q_const + w_theta_q[n_theta_q] * I_psi_const;
             } // for ( int n_theta_q = 0 ; n_theta_q < N2; n_theta_q++ )
             I_theta_q_const = J_theta_q * I_theta_q_const;
             //
             I_theta_p_const = I_theta_p_const + w_theta_p[n_theta_p] * I_theta_q_const;
         } // for ( int n_theta_p = 0 ; n_theta_p < N1 ; n_theta_p++ )


         I_subdomain_const[m-1] = J_theta_p * I_theta_p_const;
		 

     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     // end for loops
     //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

    for (int index = 0; index < 4; index++)
	{
		I += I_subdomain_const[index];
	}
     // Final output
    return I;
    

}
