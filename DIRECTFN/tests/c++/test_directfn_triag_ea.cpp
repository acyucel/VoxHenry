#include <iostream>
#include <iomanip>
#include "directfn_algorithm_ea.h"
#include "directfn_contour.h"
#include "directfn_kernel_tri.h"
#include "directfn_defs.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Triangular_EA;
using Directfn::TriangularKernel_Constant_EA;
using Directfn::TriangularKernel_RWG_WS;
using Directfn::TriangularKernel_RWG_SS;
using Directfn::TriangularKernel_nxRWG_SS;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
#ifndef M_PI
using Directfn::M_PI;
#endif
///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>  void  test_triag_ea() noexcept {

	const double  k0wn = 2 * M_PI;
	const double d = 0.1;
	const size_t  N = 10;

	double r1[3] = { 0.0, 0.0, 0.0 };
	double r2[3] = { d, 0.0, 0.0 };
	double r3[3] = { d, d, 0.0 };
	double r4[3] = { d, 0.0, d };

    SingularContour3xn cntr4p;
    cntr4p.set_points(r1, r2, r3, r4);

    unique_ptr<Triangular_EA<ParticularKernel>> up_trg_ea(new Triangular_EA<ParticularKernel>());

    up_trg_ea->set_wavenumber(k0wn);
    up_trg_ea->set(cntr4p);
    up_trg_ea->set_Gaussian_orders_4(N, N, N, N);

    up_trg_ea->calc_Iss();
    const dcomplex * ref_val = up_trg_ea->Iss();

    for (size_t k = 0; k < up_trg_ea->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}


template void test_triag_ea<TriangularKernel_Constant_EA>();
template void test_triag_ea<TriangularKernel_RWG_WS>();
template void test_triag_ea<TriangularKernel_RWG_SS>();
template void test_triag_ea<TriangularKernel_nxRWG_SS>();

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

    test_triag_ea<TriangularKernel_Constant_EA>();
    test_triag_ea<TriangularKernel_RWG_WS>();
    test_triag_ea<TriangularKernel_RWG_SS>();
    test_triag_ea<TriangularKernel_nxRWG_SS>();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file

