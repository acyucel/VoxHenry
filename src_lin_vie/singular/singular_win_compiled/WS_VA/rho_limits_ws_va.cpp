#include "constants.h"
#include "quadric_ws_va.h"



using namespace Quadrilateral_VA;

namespace Quadrilateral_VA
{
	void rho_p_limit(int argument,double theta_p, double *L_p)
	{
		switch(argument)
		{
		case(1):
			*L_p = double(2.0)/cos(theta_p);
		break;
		case(2):
			*L_p = double(2.0)/cos(theta_p);
		break;
		case(3):
			*L_p = double(2.0)/sin(theta_p);
		break;
		case(4):
			*L_p = double(2.0)/sin(theta_p);
		}
	}

	void rho_q_limit(int argument,double theta_q, double *L_q)
	{
		switch(argument)
		{
		case(1):
			*L_q = double(2.0)/cos(theta_q);
		break;
		case(2):
			*L_q = double(2.0)/sin(theta_q);
		break;
		case(3):
			*L_q = double(2.0)/cos(theta_q);
		break;
		case(4):
			*L_q = double(2.0)/sin(theta_q);

		}
	}
}