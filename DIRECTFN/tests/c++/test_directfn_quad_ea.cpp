#include <iostream>
#include <iomanip>
#include "directfn_defs.h"
#include "directfn_contour.h"
#include "directfn_algorithm_ea.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

#ifndef M_PI
using Directfn::M_PI;
#endif

using Directfn::Quadrilateral_EA;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_PlanarScalar;
using Directfn::QuadrilateralKernel_PlanarVectorWS;
using Directfn::QuadrilateralKernel_PlanarVectorSS;

using Directfn::QuadrilateralKernel_CurvilinearScalar;
using Directfn::QuadrilateralKernel_CurvilinearVectorWS;
using Directfn::QuadrilateralKernel_CurvilinearVectorSS;

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>  void  test_quad_ea_plan() noexcept {

	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N = 10;

	double r1[] = { 0.0 , 0.0 , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { d,  d , 0.0 };
	double r4[] = { 0.0 , d , 0.0 };
	double r5[] = { 0.0, 2 * d, 0.0 };
	double r6[] = { d, 2 * d, 0.0 };

    SingularContour3xn cntr6p;
    cntr6p.set_points(r1, r2, r3, r4, r5, r6);

    unique_ptr<Quadrilateral_EA<ParticularKernel>> up_quad_ea(new Quadrilateral_EA<ParticularKernel>());
    up_quad_ea->set_wavenumber(k0wn);
    up_quad_ea->set(cntr6p);
    up_quad_ea->set_Gaussian_orders_4(N, N, N, N);

    up_quad_ea->calc_Iss();
    const dcomplex * ref_val = up_quad_ea->Iss();

    for (size_t k = 0; k < up_quad_ea->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}

template <typename ParticularKernel> void test_quad_ea_curv() noexcept {

    const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N = 10;

	double dd = d / sqrt(double(2.0));

	double r1[] = { d, 0.0, 0.0 };
	double r2[] = { d, 0.0, 0.5*d };
	double r3[] = { d, 0.0, d };
	double r4[] = { dd, dd, 0.0 };
	double r5[] = { dd, dd, 0.5*d };
	double r6[] = { dd, dd, d };
	double r7[] = { 0.0, d, 0.0 };
	double r8[] = { 0.0, d, 0.5*d };
	double r9[] = { 0.0, d, d };
	double r10[] = { -dd, dd, 0.0 };
	double r11[] = { -dd, dd, 0.5*d };
	double r12[] = { -dd, dd, d };
	double r13[] = { -d, 0.0, 0.0 };
	double r14[] = { -d, 0.0, 0.5*d };
	double r15[] = { -d, 0.0, d };

    SingularContour3xn cntr6p;
    cntr6p.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15);

    unique_ptr<Quadrilateral_EA<ParticularKernel>> up_quad_ea(new Quadrilateral_EA<ParticularKernel>());
    up_quad_ea->set_wavenumber(k0wn);
    up_quad_ea->set(cntr6p);
    up_quad_ea->set_Gaussian_orders_4(N, N, N, N);

    up_quad_ea->calc_Iss();
    const dcomplex * ref_val = up_quad_ea->Iss();

    for (size_t k = 0; k < up_quad_ea->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

    cout << "plan ea scal:" << endl;
    test_quad_ea_plan<QuadrilateralKernel_PlanarScalar>();

    cout << "plan ea vect ws:" << endl;
    test_quad_ea_plan<QuadrilateralKernel_PlanarVectorWS>();

    cout << "plan ea vect ss:" << endl;
    test_quad_ea_plan<QuadrilateralKernel_PlanarVectorSS>();

    cout << "curv ea scal:" << endl;
    test_quad_ea_curv<QuadrilateralKernel_CurvilinearScalar>();

    cout << "curv ea vect ws:" << endl;
    test_quad_ea_curv<QuadrilateralKernel_CurvilinearVectorWS>();

    cout << "curv ea vect ss:" << endl;
    test_quad_ea_curv<QuadrilateralKernel_CurvilinearVectorSS>();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





