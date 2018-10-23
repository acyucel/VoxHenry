#include <iostream>
#include <iomanip>
#include "directfn_defs.h"
#include "directfn_contour.h"
#include "directfn_algorithm_st.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"

#ifndef M_PI
using Directfn::M_PI;
#endif

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Quadrilateral_ST;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_PlanarScalar;
using Directfn::QuadrilateralKernel_PlanarVectorWS;
using Directfn::QuadrilateralKernel_PlanarVectorSS;


using Directfn::QuadrilateralKernel_CurvilinearScalar;
using Directfn::QuadrilateralKernel_CurvilinearVectorWS;
using Directfn::QuadrilateralKernel_CurvilinearVectorSS;


template <typename ParticularKernel>  void  test_quad_st_plan() noexcept {

	const double k0wn = 2 * M_PI;
    const size_t  N = 10;

	const double d = 0.1;

	double r1[] = { 0.0 , 0.4*d , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { 0.8*d,  d , 0.0 };
	double r4[] = { 0.2*d , 0.8*d , 0.0 };


    SingularContour3xn cntr4p;
    cntr4p.set_points(r1, r2, r3, r4);

    unique_ptr<Quadrilateral_ST<ParticularKernel>> up_quad_st(new Quadrilateral_ST<ParticularKernel>());
    up_quad_st->set_wavenumber(k0wn);
    up_quad_st->set_Gaussian_orders_4(N, N, N, N);
    up_quad_st->set(cntr4p);

    up_quad_st->calc_Iss();
    const dcomplex * ref_val = up_quad_st->Iss();

    for (size_t k = 0; k < up_quad_st->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}

//template void test_quad_st_plan<QuadrilateralKernel_PlanarScalar>();
//template void test_quad_st_plan<QuadrilateralKernel_PlanarVectorWS>();

template <typename ParticularKernel>  void  test_quad_st_curv() noexcept {

	const double k0wn = 2.0*M_PI;
   
	const double d = 0.1;
    const size_t  N = 10;

	double dd = d / sqrt(double(2.0));

	double r1[] = { 0.0, d,0.0 };
	double r2[] = { d,0.0,0.0 };
	double r3[] = { d,0.0,d };
	double r4[] = { 0.0, d, d };
	double r5[] = { dd,dd,0.0 };
	double r6[] = { d, 0.0, 0.5*d };
	double r7[] = { dd,dd,d };
	double r8[] = { 0.0,d,0.5*d };
	double r9[] = { dd,dd,0.5*d };


    SingularContour3xn cntr9p;
    cntr9p.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9);

    unique_ptr<Quadrilateral_ST<ParticularKernel>> up_quad_st(new Quadrilateral_ST<ParticularKernel>());
    up_quad_st->set_wavenumber(k0wn);
    up_quad_st->set_Gaussian_orders_4(N, N, N, N);
    up_quad_st->set(cntr9p);

    up_quad_st->calc_Iss();
    const dcomplex * ref_val = up_quad_st->Iss();

    for (size_t k = 0; k < up_quad_st->kernel_size(); ++k) {
        cout << setprecision(17) << ref_val[k] << endl;
    }
}

//template void test_quad_st_curv<QuadrilateralKernel_CurvilinearScalar>();
//template void test_quad_st_curv<QuadrilateralKernel_CurvilinearVectorWS>();

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

    cout << "plan st scal:" << endl;
    test_quad_st_plan<QuadrilateralKernel_PlanarScalar>();

    cout << "plan st vect ws:" << endl;
    test_quad_st_plan<QuadrilateralKernel_PlanarVectorWS>();

	cout << "plan st vect ss:" << endl;
	test_quad_st_plan<QuadrilateralKernel_PlanarVectorSS>();

    cout << "curv st scal:" << endl;
    test_quad_st_curv<QuadrilateralKernel_CurvilinearScalar>();

    cout << "curv st vect ws:" << endl;
    test_quad_st_curv<QuadrilateralKernel_CurvilinearVectorWS>();

	cout << "curv st vect ss:" << endl;
	test_quad_st_curv<QuadrilateralKernel_CurvilinearVectorSS>();

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file


