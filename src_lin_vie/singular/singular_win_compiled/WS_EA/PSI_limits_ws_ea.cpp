#include "constants.h"
#include "QUADRIC_WS_EA.h"
using namespace Quadrilateral_EA;

namespace Quadrilateral_EA
{
	void PSI_limits_ws_ea ( int argument, double theta, double *psi_A, double *psi_B )
	{
    double Psi_cos = atan( cos(theta) );
	double Psi_sin = atan ( sin (theta) );

	switch(argument)
	{
                case 1: // 
                         *psi_A = double(0.0);
                         *psi_B = Psi_cos;
			 break;
                case 2: // 
                         *psi_A = Psi_cos;
                         *psi_B = M_PI / double(2.0);
			 break;
                case 3: // 
                         *psi_A = double(0.0);
                         *psi_B = Psi_sin;
			 break;
                case 4: // 
                         *psi_A = Psi_sin;
                         *psi_B = M_PI / double(2.0);
			 break;
                case 5: // 
                         *psi_A = double(0.0);
                         *psi_B = Psi_sin;
			 break;
             
				case 6: // 
                         *psi_A = Psi_sin;
                         *psi_B = M_PI / double(2.0);
			 break;
				case 7: //
						 *psi_A = double(0.0);
                         *psi_B = - Psi_cos;
			 break;
				case 8: //
						 *psi_A = - Psi_cos;
						 *psi_B = M_PI / double(2.0);
			 break;
			
		}
	}
}
