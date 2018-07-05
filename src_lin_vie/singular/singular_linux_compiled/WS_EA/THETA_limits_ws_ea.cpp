#include "constants.h"
#include "QUADRIC_WS_EA.h"
// ***********************************************************************
//			IMPLEMENTATION OF void PSI_limits
// ***********************************************************************
using namespace Quadrilateral_EA;

namespace Quadrilateral_EA
{
	void THETA_limits_ws_ea ( int argument, double *theta_A, double *theta_B )
	{
    switch(argument)
        {
                case 1: // 
                         *theta_A = double(0.0);
                         *theta_B = M_PI / double(4.0);
                         break;
                case 2: //
					     *theta_A = double(0.0);
                         *theta_B = M_PI / double(4.0);
                         break;
                case 3: // 
                         *theta_A = M_PI / double(4.0);
                         *theta_B = M_PI / double(2.0);
                         break;
                case 4: // 
                         *theta_A = M_PI / double(4.0);
                         *theta_B = M_PI / double(2.0);
                         break;
                case 5: // 
                         *theta_A = M_PI / double(2.0);
                         *theta_B = double(3.0)* M_PI / double(4.0);
                         break;
				case 6: //
                         *theta_A = M_PI / double(2.0);
                         *theta_B = double(3.0)* M_PI / double(4.0);
                         break;
				case 7: //
						 *theta_A = double(3.0)* M_PI / double(4.0);
						 *theta_B = M_PI;
						 break;
				case 8: //
						 *theta_A = double(3.0)* M_PI / double(4.0);
						 *theta_B = M_PI;
						 break;


        }
	}
}

