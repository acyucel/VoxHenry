#include "constants.h"
#include "quadric_ws_st.h"
using namespace Quadrilateral_ST;

namespace Quadrilateral_ST
{
	void lambda_limits_ws_st ( int argument, double PSI, double U, double *LAMBDA_L )
	{
		switch(argument)
		{
                case 1: // Xa
                         
                         *LAMBDA_L = ( double(1.0) - U ) / cos(PSI);
			 break;
                case 2: // Xb
                         *LAMBDA_L = ( double(1.0) - U ) / cos(PSI);
			 break;
                case 3: // Xc
                         *LAMBDA_L = double(2.0) / sin(PSI);
			 break;
                case 4: // Xd
                         *LAMBDA_L = double(2.0) / sin(PSI);
			 break;
                case 5: // Xe
                         *LAMBDA_L = ( double(1.0) + U) / cos(M_PI - PSI);
			 break;
                case 6: // Xf
                         *LAMBDA_L = ( double(1.0) + U) / cos(M_PI - PSI);
			 break;

		}
	}
}
