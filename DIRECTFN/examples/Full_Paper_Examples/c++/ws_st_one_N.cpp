#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
//#include "directfn_quad.h"
#include "directfn_defs.h"
#include "directfn_contour.h"
#include "directfn_algorithm_st.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"

using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Quadrilateral_ST;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_PlanarVectorWS;

using Directfn::relative_error;
using Directfn::max_element;

#ifndef M_PI
using Directfn::M_PI;
#endif

#ifndef DBL_EPSILON
using Directfn::DBL_EPSILON;
#endif

///////////////////////////////////////////////////////////////////////////////
void  ws_st_one_N() noexcept {

	double T;
	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 20;
	const int len = 16;

	double r1[] = { 0.0 , 0.0 , 0.0 };
	double r2[] = { d , 0.0 , 0.0 };
	double r3[] = { d,  d , 0.0 };
	double r4[] = { 0.0 , d , 0.0 };

	double max_Error_0, max_Error_1, max_Error_2, max_Error_3, max_Error_4;
	double * Error = new double[len];
	const dcomplex * I_quad;
	dcomplex * I_ref_value = new dcomplex[len];



	SingularContour3xn Q; 
	
	Q.set_points(r1, r2, r3, r4);

	unique_ptr<Quadrilateral_ST<QuadrilateralKernel_PlanarVectorWS>> up_quad_st(new Quadrilateral_ST<QuadrilateralKernel_PlanarVectorWS>());

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

	/* // Output reference value
	cout << "I_ref = " << endl;
	for (size_t i = 0; i < up_quad_st->kernel_size(); ++i) {
	cout << setprecision(20) << up_quad_st->Iss_arr(i) << endl;
	}
	cout << endl;
	*/

	std::ofstream myfile;
	myfile.open("Results_ws_st_one_N.txt");

	const int Counter = 25;

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




		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error_0 = max_element(Error, len);

		////////////////////////////N 20 20 20 //////////////////////////////////

		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(N, 20, 20, 20);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();




		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error_1 = max_element(Error,len);

		////////////////////////////20 N 20 20 //////////////////////////////////

		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(20, N, 20, 20);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();




		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error_2 = max_element(Error,len);

		////////////////////////////20 20 N 20 //////////////////////////////////

		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(20, 20, N, 20);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();




		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error_3 = max_element(Error,len);

		////////////////////////////20 20 20 N //////////////////////////////////

		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(20, 20, 20, N);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();




		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error_4 = max_element(Error,len);



		// Writing results to file
		myfile << N << " " << setprecision(20) << max_Error_0 << " " << setprecision(20) << max_Error_1 << " " << setprecision(20) << max_Error_2 << " ";
		myfile << setprecision(20) << max_Error_3 << " " << setprecision(20) << max_Error_4 << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;


}

int main(int, char *[]) {

	ws_st_one_N();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





