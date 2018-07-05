
#include <math.h>
#include "QUADRIC_WS_EA.h"



// ***********************************************************************
using namespace Quadrilateral_EA;

namespace Quadrilateral_EA
{
	void Simplex_ws_ea(double *u_p, double *v_p, double *u_q, double *v_q, double Lambda, double Psi, double Theta, double u_ )
	{
    // E_P square
    double up = u_;
    *u_p = up;
    *v_p = Lambda * sin(Psi) - 1.0;
    

    // E_Q square
    
    double uq = Lambda * cos(Psi) * cos(Theta) - u_;

  
    *u_q = uq;
    *v_q = Lambda * cos(Psi) * sin(Theta) - 1.0;
	}
}