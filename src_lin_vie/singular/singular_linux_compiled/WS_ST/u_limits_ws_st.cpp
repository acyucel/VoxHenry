#include "quadric_ws_st.h"
#include "constants.h"

using namespace Quadrilateral_ST;

namespace Quadrilateral_ST
{
	void u_limits_ws_st ( int argument, double PSI , double *U_A, double *U_B )
	{
    double u_1_psi,u_2_psi;
        switch(argument)
        {
                case 1: // Xa
                         *U_A = double(-1.0);
                         *U_B = double(1.0);
                         break;
                case 2: // Xb
					     u_1_psi = double(2.0)*tan(PSI - M_PI/double(2.0)) + double(1.0);
                         *U_A = u_1_psi;
                         *U_B = double(1.0);
                         break;
                case 3: // Xc
                         u_1_psi = double(2.0)*tan(PSI - M_PI/double(2.0)) + double(1.0);

                         *U_A = double(-1.0);
                         *U_B = u_1_psi;
                         break;
                case 4: // Xd
                         u_2_psi = double(2.0)*tan(PSI - M_PI/double(2.0)) - double(1.0);
                         *U_A = u_2_psi;
                         *U_B = double(1.0);
                         break;
                case 5: // Xe
                         u_2_psi = double(2.0)*tan(PSI - M_PI/double(2.0)) - double(1.0);
                         *U_A = double(-1.0);
                         *U_B = u_2_psi;
                         break;
                case 6: // Xf
                         *U_A = double(-1.0);
                         *U_B = double(1.0);
                         break;

        }
	}
}

