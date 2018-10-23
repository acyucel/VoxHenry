#include <iostream>
#include <iomanip>
#include "directfn_contour.h"
#include "directfn_defs.h"
#include "directfn_algorithm_va.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"

#ifndef M_PI
using Directfn::M_PI;
#endif

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Quadrilateral_VA;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;

using Directfn::QuadrilateralKernel_PlanarScalar;
using Directfn::QuadrilateralKernel_PlanarVectorWS;
using Directfn::QuadrilateralKernel_PlanarVectorSS;

using Directfn::QuadrilateralKernel_CurvilinearScalar;
using Directfn::QuadrilateralKernel_CurvilinearVectorWS;
using Directfn::QuadrilateralKernel_CurvilinearVectorSS;

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel> void test_quad_va_plan() noexcept {

	const double k0wn = 2 * M_PI;
	const size_t  N = 10;

	const double d = 0.1;

	double r1[] = { 0.0 , 0.0 , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { d,  d , 0.0 };
	double r4[] = { 0.0 , d , 0.0 };
	double r5[] = { 2 * d, d, 0.0 };
	double r6[] = { 2 * d, 2 * d, 0.0 };
	double r7[] = { d, 2 * d, 0.0 };

    SingularContour3xn cntr7p;
    cntr7p.set_points(r1, r2, r3, r4, r5, r6, r7);

    unique_ptr<Quadrilateral_VA<ParticularKernel>> up_quad_va(new Quadrilateral_VA<ParticularKernel>());
    up_quad_va->set_wavenumber(k0wn);
    up_quad_va->set(cntr7p);
    up_quad_va->set_Gaussian_orders_4(N, N, N, N);

    up_quad_va->calc_Iss();
    const dcomplex * ref_val = up_quad_va->Iss();

    for (size_t k = 0; k < up_quad_va->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}

template <typename ParticularKernel> void test_quad_va_curv() noexcept {

	const double k0wn = 2 * M_PI;
    const size_t  N = 10;
	const double d = 0.1;

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
	double r10[] = { 0.0, d, 1.5*d };
	double r11[] = { 0.0, d, 2.0*d };
	double r12[] = { -dd, dd, d };
	double r13[] = { -dd, dd, 1.5*d };
	double r14[] = { -dd, dd, 2.0*d };
	double r15[] = { -d, 0.0, d };
	double r16[] = { -d, 0.0, 1.5*d };
	double r17[] = { -d, 0.0, 2.0*d };

    SingularContour3xn cntr7p;
    cntr7p.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17);

    unique_ptr<Quadrilateral_VA<ParticularKernel>> up_quad_va(new Quadrilateral_VA<ParticularKernel>());
    up_quad_va->set_wavenumber(k0wn);
    up_quad_va->set(cntr7p);
    up_quad_va->set_Gaussian_orders_4(N, N, N, N);

    up_quad_va->calc_Iss();
    const dcomplex * ref_val = up_quad_va->Iss();

    for (size_t k = 0; k < up_quad_va->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

    cout << "plan va scal:" << endl;
    test_quad_va_plan<QuadrilateralKernel_PlanarScalar>();

    cout << "plan va vect ws:" << endl;
    test_quad_va_plan<QuadrilateralKernel_PlanarVectorWS>();

    cout << "plan va vect ss:" << endl;
    test_quad_va_plan<QuadrilateralKernel_PlanarVectorSS>();

    cout << "curv va scal:" << endl;
    test_quad_va_curv<QuadrilateralKernel_CurvilinearScalar>();

    cout << "curv va vect ws:" << endl;
    test_quad_va_curv<QuadrilateralKernel_CurvilinearVectorWS>();

    cout << "curv va vect ss:" << endl;
    test_quad_va_curv<QuadrilateralKernel_CurvilinearVectorSS>();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





