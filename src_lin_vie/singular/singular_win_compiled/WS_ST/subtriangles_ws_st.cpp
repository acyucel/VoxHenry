#include "quadric_ws_st.h"

using namespace Quadrilateral_ST;

namespace Quadrilateral_ST
{
	void subtriangles_ws_st(double U_sub, double V_sub, int argument, double *U, double *V)
	{
		switch(argument)
		{
		case 1:
			*U  = U_sub;
			*V  = V_sub;
			break;
		case 2:
			*U = - V_sub;
			*V  =  U_sub;
			break;
		case 3:
			*U = - U_sub;
			*V = - V_sub;
			break;
		case 4:
			*U  = V_sub;
			*V = - U_sub;
			break;
		}
	}
}
