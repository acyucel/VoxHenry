#include <iostream>
#include <iomanip>
#include "directfn_algorithm_st.h"
#include "directfn_contour.h"
#include "directfn_kernel_tri.h"
#include "directfn_defs.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Triangular_ST;
using Directfn::TriangularKernel_Constant_ST;
using Directfn::TriangularKernel_RWG_WS;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
#ifndef M_PI
using Directfn::M_PI;
#endif

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>  void  test_triag_st() noexcept {

	const double  k0wn = 2 * M_PI;
	const double d = 0.1;
	const size_t  N = 10;

	double r1[3] = { 0.0, 0.0, 0.0 };
	double r2[3] = { d, 0.0, 0.0 };
	double r3[3] = { 0.0, d, 0.0 };

    SingularContour3xn cntr3p;
    cntr3p.set_points(r1, r2, r3);

    unique_ptr<Triangular_ST<ParticularKernel>> up_trg_st(new Triangular_ST<ParticularKernel>());

    up_trg_st->set_wavenumber(k0wn);
    up_trg_st->set(cntr3p);
    up_trg_st->set_Gaussian_orders_4(N, N, N, N);

    up_trg_st->calc_Iss();
    const dcomplex * ref_val = up_trg_st->Iss();

    for (size_t k = 0; k < up_trg_st->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}


template void test_triag_st<TriangularKernel_Constant_ST>();
template void test_triag_st<TriangularKernel_RWG_WS>();

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

    test_triag_st<TriangularKernel_Constant_ST>();
    test_triag_st<TriangularKernel_RWG_WS>();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





