#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
//#include "directfn_quad.h'
//#include "directfn_tri.h"
#include "directfn_contour.h"
#include "directfn_defs.h"
#include "directfn_algorithm_st.h"
#include "directfn_kernel_quad_scal.h"
#include <string>

using  std::cout;
using  std::endl;
using  std::setprecision;
using  std::string;

using Directfn::Quadrilateral_ST;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_PlanarScalar;


#ifndef M_PI
using Directfn::M_PI;
#endif

#ifndef DBL_EPSILON
using Directfn::DBL_EPSILON;
#endif
///////////////////////////////////////////////////////////////////////////////
 void  shape_test_st(int m, double step) noexcept {

    const double k0wn = 2 * M_PI;
    
    const double d = 0.1;
	const size_t  N_ref = 30;

    double r1[] = { 0.0 , 0.0 , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { d + (m-1)*step*d,  d , 0.0 };
	double r4[] = { 0.0 + (m-1)*step*d, d , 0.0 };

	double Error_quad;
	const dcomplex * I_quad;
	dcomplex value_quad;



    SingularContour3xn Q; // Quadrilateral
    Q.set_points(r1, r2, r3, r4);

    unique_ptr<Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>> up_quad_st(new Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>());

	// Setting parameters
	up_quad_st->set_wavenumber(k0wn);
	up_quad_st->set(Q);
	up_quad_st->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values for ST case ..." << endl;
	up_quad_st->calc_Iss(); 
	cout << "Reference values computed "  << endl;

	const dcomplex * I_ref = up_quad_st->Iss();
	dcomplex ref_val = I_ref[0];


	cout << setprecision(17) << ref_val << endl;
	

	std::ofstream myfile;
	char* filename;
	sprintf(filename, "Results_different_shape_%d.txt", m);
	myfile.open(filename);

	const int Counter = 30;

	cout << "Convergence test for ST case starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(N, N, N, N);

		// DIRECTFN-quad
		//std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
		up_quad_st->calc_Iss();
		//std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

		I_quad = up_quad_st->Iss();
		value_quad = I_quad[0];
		Error_quad = fabs(abs((value_quad - ref_val)) / abs(ref_val) + DBL_EPSILON);
		//T_quad = std::chrono::duration_cast<std::chrono::microseconds>(end_t_quad - begin_t_quad).count();

		cout << "N: " << N << endl;
		//cout << "Runtime_quad: " << setprecision(4) << T_quad << " [mcsec]" << endl;
		cout << "I_ST_quad = " << setprecision(20) << value_quad << endl;


		// Writing results to file
		myfile << N << " ";
		myfile << setprecision(20) << Error_quad << " " << endl;
	}

	myfile.close();
	cout << "Convergence test for ST case completed." << endl;

}

///////////////////////////////////////////////////////////////////////////////

int main(int , char *  []) {

	double step = 0.229;
	for (int m = 1; m <= 8; m++)
	{
		shape_test_st(m, step);
	}
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





