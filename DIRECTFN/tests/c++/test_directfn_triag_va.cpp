#include <iostream>
#include <iomanip>
#include "directfn_algorithm_va.h"
#include "directfn_contour.h"
#include "directfn_kernel_tri.h"
#include "directfn_defs.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Triangular_VA;
using Directfn::TriangularKernel_Constant_VA;
using Directfn::TriangularKernel_RWG_WS;
using Directfn::TriangularKernel_RWG_SS;
using Directfn::TriangularKernel_nxRWG_SS;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
#ifndef M_PI
using Directfn::M_PI;
#endif

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>  void  test_triag_va() noexcept {

	const double  k0wn = 2 * M_PI;
	const double d = 0.1;
	const size_t  N = 10;

	double r1[3] = { 0.0, d, 0.0 };
	double r2[3] = { 0.0, 0.0, 0.0 };
	double r3[3] = { d, 0.0, 0.0 };
	double r4[3] = { d, d, d };
	double r5[3] = { 2*d, d, d };

    SingularContour3xn cntr5p;
    cntr5p.set_points(r1, r2, r3, r4, r5);

    unique_ptr<Triangular_VA<ParticularKernel>> up_trg_va(new Triangular_VA<ParticularKernel>());

    up_trg_va->set_wavenumber(k0wn);
    up_trg_va->set(cntr5p);
    up_trg_va->set_Gaussian_orders_4(N, N, N, N);

    up_trg_va->calc_Iss();
    const dcomplex * ref_val = up_trg_va->Iss();

    for (size_t k = 0; k < up_trg_va->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}


template void test_triag_va<TriangularKernel_Constant_VA>();
template void test_triag_va<TriangularKernel_RWG_WS>();
template void test_triag_va<TriangularKernel_RWG_SS>();
template void test_triag_va<TriangularKernel_nxRWG_SS>();

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

    test_triag_va<TriangularKernel_Constant_VA>();
    test_triag_va<TriangularKernel_RWG_WS>();
    test_triag_va<TriangularKernel_RWG_SS>();
    test_triag_va<TriangularKernel_nxRWG_SS>();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file






