#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
//#include "directfn_quad.h"
#include "directfn_defs.h"
#include "directfn_contour.h"
#include "directfn_algorithm_st.h"
#include "directfn_kernel_quad_scal.h"


using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Quadrilateral_ST;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_CurvilinearScalar;

using Directfn::relative_error;
using Directfn::max_element;

#ifndef M_PI
using Directfn::M_PI;
#endif

#ifndef DBL_EPSILON
using Directfn::DBL_EPSILON;
#endif

///////////////////////////////////////////////////////////////////////////////
void  st_curv() noexcept {

	double T;
	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 20;
	const int len = 1;

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

	double Error;
	const dcomplex * I_quad;
	dcomplex * I_ref_value = new dcomplex[len];



	SingularContour3xn Q;

	Q.set_points(r1, r5, r2, r8, r9, r6, r4, r7, r3);

	unique_ptr<Quadrilateral_ST<QuadrilateralKernel_CurvilinearScalar>> up_quad_st(new Quadrilateral_ST<QuadrilateralKernel_CurvilinearScalar>());

	// Setting parameters
	up_quad_st->set_wavenumber(k0wn);
	up_quad_st->set(Q);
	up_quad_st->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values ..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_st->calc_Iss();
	std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference value computed in " << T << " milliseconds" << endl;

	up_quad_st->copy_Iss_array_values_to(I_ref_value);

	 // Output reference value
	cout << "I_ref = ";
	for (size_t i = 0; i < up_quad_st->kernel_size(); ++i) {
	cout << setprecision(20) << up_quad_st->Iss_arr(i) << endl;
	}
	cout << endl;
	

	std::ofstream myfile;
	myfile.open("Results_st_curv.txt");

	const int Counter = 25;

	cout << "Convergence test starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{

		cout << "N: " << N << endl;
		cout << endl;

		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(N, N, N, N);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();




		Error = relative_error(I_quad[0], I_ref_value[0]);


		// Writing results to file
		myfile << N << " " << setprecision(20) << Error << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;


}

int main(int, char *[]) {

	st_curv();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





