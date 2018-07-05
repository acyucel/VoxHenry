#include "constants.h"
#include "quadric_ws_va.h"

using namespace Quadrilateral_VA;

namespace Quadrilateral_VA
{
	void theta_q_limits ( int argument, double *theta_A, double *theta_B )
	{
		switch(argument)
		{
                case 1: // 
                         *theta_A = double(0.0);
                         *theta_B = M_PI / double(4.0);
			 break;
                case 2: // 
                        *theta_A = M_PI / double(4.0);
                         *theta_B = M_PI / double(2.0);
			 break;
                case 3: // 
                         *theta_A = double(0.0);
                         *theta_B = M_PI / double(4.0);
			 break;
                case 4: // 
                         *theta_A = M_PI / double(4.0);
                         *theta_B = M_PI / double(2.0);
			 break;

		}
	}
}
