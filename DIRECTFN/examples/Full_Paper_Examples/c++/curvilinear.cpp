#include <iostream>
#include <iomanip>
#include <chrono>
#include <fstream>
#include <cfloat>
//#include "directfn_quad.h"
#include "directfn_defs.h"
#include "directfn_contour.h"
#include "directfn_algorithm_st.h"
#include "directfn_algorithm_ea.h"
#include "directfn_algorithm_va.h"
#include "directfn_kernel_quad_vect.h"


using  std::cout;
using  std::endl;
using  std::setprecision;

using Directfn::Quadrilateral_ST;
using Directfn::Quadrilateral_EA;
using Directfn::Quadrilateral_VA;
using Directfn::SingularContour3xn;
using Directfn::dcomplex;
using Directfn::QuadrilateralKernel_CurvilinearVectorWS;
using Directfn::QuadrilateralKernel_CurvilinearVectorSS;

using Directfn::relative_error;
using Directfn::max_element;

#ifndef M_PI
using Directfn::M_PI;
#endif

//#ifndef DBL_EPSILON
//using Directfn::DBL_EPSILON;
//#endif

///////////////////////////////////////////////////////////////////////////////
void  st_curv() noexcept {

	double T;
	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 20;
	const int len = 16;

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

	double * Error = new double[len];
	double max_Error;
	const dcomplex * I_quad;
	dcomplex * I_ref_value = new dcomplex[len];



	SingularContour3xn Q;

	Q.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9);

	unique_ptr<Quadrilateral_ST<QuadrilateralKernel_CurvilinearVectorWS>> up_quad_st(new Quadrilateral_ST<QuadrilateralKernel_CurvilinearVectorWS>());

	// Setting parameters
	up_quad_st->set_wavenumber(k0wn);
	up_quad_st->set(Q);
	up_quad_st->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values for ST case..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_st->calc_Iss();
	std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference values computed in " << T << " milliseconds" << endl;

	up_quad_st->copy_Iss_array_values_to(I_ref_value);

	 // Output reference value
	cout << "I_ref = ";
	for (size_t i = 0; i < up_quad_st->kernel_size(); ++i) {
	cout << setprecision(17) << up_quad_st->Iss_arr(i) << endl;
	}
	cout << endl;
	

	std::ofstream myfile;
	myfile.open("Results_st_curv.txt");

	const int Counter = 25;

	cout << "Convergence test starting for ST case..." << endl;
	for (int N = 1; N <= Counter; N++)
	{

		cout << "N: " << N << endl;
		cout << endl;

		// Setting Gaussian orders
		up_quad_st->set_Gaussian_orders_4(N, N, N, N);
		up_quad_st->calc_Iss();
		I_quad = up_quad_st->Iss();


		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error = max_element(Error, len);



		// Writing results to file
		myfile << N << " " << setprecision(17) << max_Error << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;


}


// Edge adjacent case
void  ea_curv() noexcept {

	double T;
	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 20;
	const int len = 16;

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

	double * Error = new double[len];
	double max_Error;
	const dcomplex * I_quad;
	dcomplex * I_ref_value = new dcomplex[len];



	SingularContour3xn Q;

	Q.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15);

	unique_ptr<Quadrilateral_EA<QuadrilateralKernel_CurvilinearVectorSS>> up_quad_ea(new Quadrilateral_EA<QuadrilateralKernel_CurvilinearVectorSS>());

	// Setting parameters
	up_quad_ea->set_wavenumber(k0wn);
	up_quad_ea->set(Q);
	up_quad_ea->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values for EA case..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_ea->calc_Iss();
	std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference value computed in " << T << " milliseconds" << endl;

	up_quad_ea->copy_Iss_array_values_to(I_ref_value);

	// Output reference value
	cout << "I_ref = ";
	for (size_t i = 0; i < up_quad_ea->kernel_size(); ++i) {
		cout << setprecision(17) << up_quad_ea->Iss_arr(i) << endl;
	}
	cout << endl;


	std::ofstream myfile;
	myfile.open("Results_ea_curv.txt");

	const int Counter = 25;

	cout << "Convergence test starting for EA case ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{

		cout << "N: " << N << endl;
		cout << endl;

		// Setting Gaussian orders
		up_quad_ea->set_Gaussian_orders_4(N, N, N, N);
		up_quad_ea->calc_Iss();
		I_quad = up_quad_ea->Iss();

		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error = max_element(Error, len);


		// Writing results to file
		myfile << N << " " << setprecision(17) << max_Error << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;


}

// Vertex adjacent case
void  va_curv() noexcept {

	double T;
	const double k0wn = 2 * M_PI;

	const double d = 0.1;
	const size_t  N_ref = 20;
	const int len = 16;

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
	double r16[] = { 0.0, d, 1.5*d };
	double r17[] = { 0.0, d, 2.0*d };
	double r18[] = { -dd, dd, d };
	double r19[] = { -dd, dd, 1.5*d };
	double r20[] = { -dd, dd, 2.0*d };
	double r21[] = { -d, 0.0, d };
	double r22[] = { -d, 0.0, 1.5*d };
	double r23[] = { -d, 0.0, 2.0*d };


	double * Error = new double[len];
	double max_Error;
	const dcomplex * I_quad;
	dcomplex * I_ref_value = new dcomplex[len];



	SingularContour3xn Q;

	Q.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9, r16, r17, r18, r19, r20, r21, r22, r23);

	unique_ptr<Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorSS>> up_quad_va(new Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorSS>());

	// Setting parameters
	up_quad_va->set_wavenumber(k0wn);
	up_quad_va->set(Q);
	up_quad_va->set_Gaussian_orders_4(N_ref, N_ref, N_ref, N_ref);

	// Calculating reference value
	cout << "Computing reference values for VA case..." << endl;
	std::chrono::steady_clock::time_point  begin_t_quad = std::chrono::steady_clock::now();
	up_quad_va->calc_Iss();
	std::chrono::steady_clock::time_point  end_t_quad = std::chrono::steady_clock::now();

	T = std::chrono::duration_cast<std::chrono::milliseconds>(end_t_quad - begin_t_quad).count();
	cout << "Reference value computed in " << T << " milliseconds" << endl;

	up_quad_va->copy_Iss_array_values_to(I_ref_value);

	// Output reference value
	cout << "I_ref = ";
	for (size_t i = 0; i < up_quad_va->kernel_size(); ++i) {
		cout << setprecision(17) << up_quad_va->Iss_arr(i) << endl;
	}
	cout << endl;


	std::ofstream myfile;
	myfile.open("Results_va_curv.txt");

	const int Counter = 25;

	cout << "Convergence test starting for VA case ..." << endl;
	for (int N = 1; N <= Counter; N++)
	{

		cout << "N: " << N << endl;
		cout << endl;

		// Setting Gaussian orders
		up_quad_va->set_Gaussian_orders_4(N, N, N, N);
		up_quad_va->calc_Iss();
		I_quad = up_quad_va->Iss();

		for (int i = 0; i < len; i++)
		{
			Error[i] = relative_error(I_quad[i], I_ref_value[i]);
		}

		max_Error = max_element(Error, len);


		// Writing results to file
		myfile << N << " " << setprecision(17) << max_Error << endl;
	}

	myfile.close();
	cout << "Convergence test completed." << endl;


}

int main(int, char *[]) {

	st_curv();
	ea_curv();
	va_curv();

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

// End of the file





