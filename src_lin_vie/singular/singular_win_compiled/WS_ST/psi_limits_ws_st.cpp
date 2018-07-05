#include "constants.h"
#include "quadric_ws_st.h"
using namespace Quadrilateral_ST;

namespace Quadrilateral_ST{
	void psi_limits_ws_st ( int argument, double *psi_A, double *psi_B )
	{
		switch(argument)
		{
                case 1: // Xa
                         *psi_A = double(0.0);
                         *psi_B = M_PI / double(4.0);
			 break;
                case 2: // Xb
                         *psi_A = M_PI / double(4.0);
                         *psi_B = M_PI / double(2.0);
			 break;
                case 3: // Xc
                         *psi_A = M_PI / double(4.0);
                         *psi_B = M_PI / double(2.0);
			 break;
                case 4: // Xd
                         *psi_A = M_PI / double(2.0);
                         *psi_B = double(3.0) * M_PI / double(4.0);
			 break;
                case 5: // Xe
                         *psi_A = M_PI / double(2.0);
                         *psi_B = double(3.0) * M_PI / double(4.0);
			 break;
                case 6: // Xf
                         *psi_A = double(3.0) * M_PI / double(4.0);
                         *psi_B = M_PI;
			 break;

		}
	}
}
