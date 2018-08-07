#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include "directfn_defs.h"
#include "directfn_contour.h"
#include "directfn_algorithm_st.h"
#include "directfn_kernel_quad_scal.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Quadrilateral_ST;
using Directfn::SingularContour3xn;
using Directfn::QuadrilateralKernel_PlanarScalar;

using Directfn::dcomplex;
using Directfn::relative_error;
using Directfn::max_element;

#ifndef M_PI
using Directfn::M_PI;
#endif

#ifndef DBL_EPSILON
using Directfn::DBL_EPSILON;
#endif

///////////////////////////////////////////////////////////////////////////////
void  st_one_param_const() noexcept {

	double T;

	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 20;
	const int len = 1;

	const double r1[] = { 0.0 , 0.4*d , 0.0 };
	const double r2[] = { d , 0.0 , 0.0 };
	const double r3[] = { 0.8*d,  d , 0.0 };
	const double r4[] = { 0.2*d , 0.8*d , 0.0 };

	double Error_0, Error_1, Error_2, Error_3, Error_4;
	const dcomplex * I_quad;
	dcomplex * I_ref_value = new dcomplex[len];



	SingularContour3xn Q;

	Q.set_points(r1, r2, r3, r4);

	unique_ptr<Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>> up_quad_st(new Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>());

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


	std::ofstream myfile;
	myfile.open("Results_st_one_param_const.txt");

	const int Counter = 20;

	cout << "Convergence test starting ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{

		cout << "N: " << N << endl;
		cout << endl;
		////////////////////////////Same N //////////////////////////////////
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(N, N, N, N);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();




		Error_0 = relative_error(I_quad[0], I_ref_value[0]);

		////////////////////////////N 20 20 20 //////////////////////////////////
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(N, 20, 20, 20);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();


		Error_1 = relative_error(I_quad[0], I_ref_value[0]);

		////////////////////////////20 N 20 20 //////////////////////////////////
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(20, N, 20, 20);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();


        Error_2 = relative_error(I_quad[0], I_ref_value[0]);
	
		////////////////////////////20 20 N 20 //////////////////////////////////
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(20, 20, N, 20);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();


		Error_3 = relative_error(I_quad[0], I_ref_value[0]);


		////////////////////////////20 20 20 N //////////////////////////////////
		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(20, 20, 20, N);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();




	    Error_4 = relative_error(I_quad[0], I_ref_value[0]);
	



		// Writing results to file
		myfile << N << " " << setprecision(20) << Error_0 << " " << setprecision(20) << Error_1 << " " << setprecision(20) << Error_2 << " ";
		myfile << setprecision(20) << Error_3 << " " << setprecision(20) << Error_4 << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;


}

int main(int, char *[]) {

	st_one_param_const();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





