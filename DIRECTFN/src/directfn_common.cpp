#include <algorithm>
#include <iostream>
#include <iomanip>
#include "directfn_common.h"
#include "directfn_interface.h"

using std::max;
using std::fabs;
using std::cout;
using std::endl;
using std::setprecision;
using std::setw;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

// Template definition
template <typename T_first_input, typename T_secnd_input, typename T_result>
T_result  vector_dot(const T_first_input x[3], const T_secnd_input y[3]) {
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
}

double vector_dot(const double x[3], const double y[3]) {
    return vector_dot<double, double, double>(x, y);
}

// Template instantiation
template double   vector_dot(const double   [3], const double   [3]);
template dcomplex vector_dot(const double   [3], const dcomplex [3]);
template dcomplex vector_dot(const dcomplex [3], const double   [3]);

///----------------------------------------------------------------------------

double sqrt_vector_dot(const double x[3], const double y[3]) {

    const double max_abs_0 = max(fabs(x[0]), fabs(y[0]));
    const double max_abs_1 = max(fabs(x[1]), fabs(y[1]));
    const double max_abs_2 = max(fabs(x[2]), fabs(y[2]));
    const double xy_ma = max(max(max_abs_0, max_abs_1), max_abs_2);

    return xy_ma * sqrt((x[0] / xy_ma) * (y[0] / xy_ma) +
                        (x[1] / xy_ma) * (y[1] / xy_ma) +
                        (x[2] / xy_ma) * (y[2] / xy_ma));
}

///----------------------------------------------------------------------------

// Template definition
template <typename Real_T, typename Complex_T>
void vector_cross(const Complex_T x[3], const Real_T y[3], Complex_T z[3]) {

    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
}

// Template instantiation
template void vector_cross(const double   [3], const double [3], double   [3]);
template void vector_cross(const dcomplex [3], const double [3], dcomplex [3]);

///----------------------------------------------------------------------------

void calc_normal(double nr[3], const double r1[3], const double r2[3], const double r3[3]) noexcept {

    double l1[3];
    double l2[3];

    for (int i = 0; i < 3; ++i) {
        l1[i] = r2[i] - r1[i];
        l2[i] = r3[i] - r1[i];
    }
    vector_cross(l1, l2, nr);

    const double n12m1 = 1.0 / sqrt(vector_dot<double, double, double>(nr, nr));
    for(int i = 0; i < 3; ++i) {
        nr[i] *= n12m1;
    }
}

//void  print(const DirectfnInterface * figure_xx,
//            const dcomplex * ref_val, const dcomplex * eth_val) noexcept {

//    cout << endl << "Test for " << figure_xx->name() << endl;
//    for (size_t i = 0; i < figure_xx->kernels_size(); ++i) {

//        cout << setprecision(17) << " eth = " << eth_val[i] << endl;
//        cout << setprecision(17) << " ref = " << ref_val[i] << endl;
//        const dcomplex err = eth_val[i] - ref_val[i];
//        const double erval = abs(err) / abs(eth_val[i]);

//        const double ERR_THRESHOLD = 1e-14;
//        if (erval > ERR_THRESHOLD) {
//            cout << "WARNING!  Too high error" << endl;
//        }
//        cout << "     err = " << erval << endl;
//    }
//}

void calc_rotation_matrix_co(const double u[3], const double theta, double r[3]) noexcept {

    const size_t nsz3 = 3;
    const size_t x = 0;
    const size_t y = 1;
    const size_t z = 2;

    double R[9];
    // Define rotation matrix
    R[0 * nsz3 + 0] = cos(theta) + u[x] * u[x] * (1.0 - cos(theta));
    R[0 * nsz3 + 1] = u[x] * u[y] * (1.0 - cos(theta)) - u[z] * sin(theta);
    R[0 * nsz3 + 2] = u[x] * u[z] * (1.0 - cos(theta)) + u[y] * sin(theta);

    R[1 * nsz3 + 0] = u[y] * u[x] * (1.0 - cos(theta)) + u[z] * sin(theta);
    R[1 * nsz3 + 1] = cos(theta) + u[y] * u[y] * (1.0 - cos(theta));
    R[1 * nsz3 + 2] = u[y] * u[z] * (1.0 - cos(theta)) - u[x] * sin(theta);

    R[2 * nsz3 + 0] = u[z] * u[x] * (1.0 - cos(theta)) - u[y] * sin(theta);
    R[2 * nsz3 + 1] = u[z] * u[y] * (1.0 - cos(theta)) + u[x] * sin(theta);
    R[2 * nsz3 + 2] = cos(theta) + u[z] * u[z] * (1.0 - cos(theta));

    // Rotate in temporary vector
    double rez[3] = {0, 0, 0};
    for (size_t i = 0; i < nsz3; ++i) {
        double sum_i = 0;
        for (size_t k = 0; k < nsz3; ++k) {
            sum_i += R[i * nsz3 + k] * r[k];
        }
        rez[i] = sum_i;
    }

    // Assign the result to the output vector
    r[0] = rez[0];
    r[1] = rez[1];
    r[2] = rez[2];
}

//double max_element(double Error[]) {

//    const size_t size_Error = sizeof(Error);
//    double max = Error[0];
//    for (size_t i = 1; i < size_Error; ++i) {
//        if (Error[i] > max) {
//            max = Error[i];
//        }
//    }
//    return max;
//}

double max_element(const double * const Error, const size_t err_len) {

    const size_t size_Error = err_len;
    double max = Error[0];
    for (size_t i = 1; i < size_Error; ++i) {
        if (Error[i] > max) {
            max = Error[i];
        }
    }
    return max;
}

double relative_error(dcomplex x, dcomplex x_ref) {

    //if ((abs(x) <= DBL_EPSILON) && (abs(x_ref)) <= DBL_EPSILON)
    //if (abs(x_ref) < DBL_EPSILON) {
        //return DBL_EPSILON;
    //}
    //else {
    //return fabs(fabs(abs(x - x_ref) / abs(x_ref)) + DBL_EPSILON);

	//if ((abs(x) <= DBL_EPSILON) && (abs(x_ref)) <= DBL_EPSILON) {
	if (abs(x_ref) <= DBL_EPSILON) {
	    return DBL_EPSILON;
	} else {
		return abs(x - x_ref) / abs(x_ref) + DBL_EPSILON;
	}
	

//	const double diff_xxref = abs(x - x_ref);
//	const double max_x_xref = std::max<double>(abs(x), abs(x_ref));
//	const double abs_x_ref = abs(x_ref);

//	if (diff_xxref < DBL_EPSILON * abs_x_ref) {
	//	return DBL_EPSILON;
//	}
//	return diff_xxref / abs_x_ref;

}

///////////////////////////////////////////////////////////////////////////////

template <>
void triag_st_subtriangles<1>(const double ETA_sub, const double XI_sub, double & ETA, double & XI) noexcept {

    ETA = ETA_sub;
    XI  = XI_sub;
}

template <>
void triag_st_subtriangles<2>(const double ETA_sub, const double XI_sub, double & ETA, double & XI) noexcept {

    const double dsqrt3 = sqrt(double(3.0));

    ETA = 0.5 * (1.0 -  ETA_sub - dsqrt3 * XI_sub);
    XI  = 0.5 * (dsqrt3  +  dsqrt3 * ETA_sub - XI_sub);
}

template <>
void triag_st_subtriangles<3>(const double ETA_sub, const double XI_sub, double & ETA, double & XI) noexcept {

    const double dsqrt3 = sqrt(double(3.0));

    ETA = 0.5 * (-1.0 - ETA_sub + dsqrt3 * XI_sub);
    XI  = 0.5 * (dsqrt3 -  dsqrt3 * ETA_sub - XI_sub);
}

///----------------------------------------------------------------------------

template <>
void  quadr_st_subtriangles<1>(const double U_sub, const double V_sub, double & Uout, double & Vout) noexcept {

    Uout = U_sub;
    Vout = V_sub;
}

template <>
void  quadr_st_subtriangles<2>(const double U_sub, const double V_sub, double & Uout, double & Vout) noexcept {

    Uout = -V_sub;
    Vout =  U_sub;
}

template <>
void  quadr_st_subtriangles<3>(const double U_sub, const double V_sub, double & Uout, double & Vout) noexcept {

    Uout = -U_sub;
    Vout = -V_sub;
}

template <>
void  quadr_st_subtriangles<4>(const double U_sub, const double V_sub, double & Uout, double & Vout) noexcept {

    Uout =  V_sub;
    Vout = -U_sub;
}

///----------------------------------------------------------------------------
/// Psi limits Triang

template <>
void triag_st_psi_limits<1>(double & psi_A, double & psi_B) noexcept {
    psi_A = double(0.0);
    psi_B = M_PI / double(3.0);
}

template <>
void triag_st_psi_limits<2>(double & psi_A, double & psi_B) noexcept {
    psi_A = M_PI / double(3.0);
    psi_B = double(2.0) * M_PI / double(3.0);
}

template <>
void triag_st_psi_limits<3>(double & psi_A, double & psi_B) noexcept {
    psi_A = double(2.0) * M_PI / double(3.0);
    psi_B = M_PI;
}

template <>
void triag_st_psi_limits<4>(double & psi_A, double & psi_B) noexcept {
    psi_A = double(2.0) * M_PI / double(3.0);
    psi_B = M_PI;
}

template <>
void triag_st_psi_limits<5>(double & psi_A, double & psi_B) noexcept {
    psi_A = double(0.0);
    psi_B = M_PI / double(3.0);
}

template <>
void triag_st_psi_limits<6>(double & psi_A, double & psi_B) noexcept {
    psi_A = double(0.0);
    psi_B = M_PI / double(3.0);
}

template <>
void triag_st_psi_limits<7>(double & psi_A, double & psi_B) noexcept {
    psi_A = M_PI / double(3.0);
    psi_B = double(2.0) * M_PI / double(3.0);
}

template <>
void triag_st_psi_limits<8>(double & psi_A, double & psi_B) noexcept {
    psi_A = double(2.0) * M_PI / double(3.0);
    psi_B = M_PI;
}

///----------------------------------------------------------------------------

/// Psi limits Quads

template <>
void quad_st_psi_limits<1>(double & psi_A, double & psi_B) noexcept {
    psi_A = 0.0;
    psi_B = 0.25 * M_PI;
}

template <>
void quad_st_psi_limits<2>(double & psi_A, double & psi_B) noexcept {
    psi_A = 0.25 * M_PI;
    psi_B = 0.5  * M_PI;
}

template <>
void quad_st_psi_limits<3>(double & psi_A, double & psi_B) noexcept {
    psi_A = 0.25 * M_PI;
    psi_B = 0.5  * M_PI;
}

template <>
void quad_st_psi_limits<4>(double & psi_A, double & psi_B) noexcept {
    psi_A = 0.5  * M_PI;
    psi_B = 0.75 * M_PI;
}

template <>
void quad_st_psi_limits<5>(double & psi_A, double & psi_B) noexcept {
    psi_A = 0.5  * M_PI;
    psi_B = 0.75 * M_PI;
}

template <>
void quad_st_psi_limits<6>(double & psi_A, double & psi_B) noexcept {
    psi_A = 0.75  * M_PI;
    psi_B = M_PI;
}

///----------------------------------------------------------------------------

template <>
void  triag_st_eta_limits<1>(const double , double & ETA_A, double & ETA_B) noexcept {
    ETA_A = double(0.0);
    ETA_B = double(1.0);
}

template <>
void  triag_st_eta_limits<2>(const double , double & ETA_A, double & ETA_B) noexcept {
    ETA_A = double(0.0);
    ETA_B = double(1.0);
}

template <>
void  triag_st_eta_limits<3>(const double PSI, double & ETA_A, double & ETA_B) noexcept {

    const double b_1 = tan(M_PI - PSI) / sqrt(double(3.0));
    ETA_A = (double(1.0) - b_1) / (double(1.0) + b_1);
    ETA_B = double(1.0);
}

template <>
void  triag_st_eta_limits<4>(const double PSI, double & ETA_A, double & ETA_B) noexcept {

    const double b_1 = tan(M_PI - PSI) / sqrt( double(3.0));
    ETA_B = ( double(1.0) - b_1 ) / ( double(1.0) + b_1);
    ETA_A = double(0.0);
}

template <>
void  triag_st_eta_limits<5>(const double PSI, double & ETA_A, double & ETA_B) noexcept {

    const double b_2   = tan(PSI) / sqrt( double(3.0) );
    ETA_A = -( double(1.0) - b_2) / (double(1.0) + b_2);
    ETA_B = double(0.0);
}

template <>
void  triag_st_eta_limits<6>(const double PSI, double & ETA_A, double & ETA_B) noexcept {

    const double b_2   = tan(PSI) / sqrt( double(3.0) );
    ETA_B = -( double(1.0) - b_2 ) / ( double(1.0) + b_2);
    ETA_A = -double(1.0);
}

template <>
void  triag_st_eta_limits<7>(const double , double & ETA_A, double & ETA_B) noexcept {
    ETA_A = -double(1.0);
    ETA_B =  double(0.0);
}

template <>
void  triag_st_eta_limits<8>(const double , double & ETA_A, double & ETA_B) noexcept {
    ETA_A = -double(1.0);
    ETA_B =  double(0.0);
}

///----------------------------------------------------------------------------

template <>
void quad_st_eta_limits<16>(const double, double & U_a, double & U_b) noexcept {

    U_a = -1.0;
    U_b =  1.0;
}

template <>
void quad_st_eta_limits<2>(const double Psi, double & U_a, double & U_b) noexcept {

    U_a = 2.0 * tan(Psi - 0.5 * M_PI) + 1.0;
    U_b = 1.0;
}

template <>
void quad_st_eta_limits<3>(const double Psi, double & U_a, double & U_b) noexcept {

    U_a = -1.0;
    U_b = 2.0 * tan(Psi - 0.5 * M_PI) + 1.0;
}

template <>
void quad_st_eta_limits<4>(const double Psi, double & U_a, double & U_b) noexcept {

    U_a = 2.0 * tan(Psi - 0.5 * M_PI) - 1.0;
    U_b = 1.0;
}

template <>
void quad_st_eta_limits<5>(const double Psi, double & U_a, double & U_b) noexcept {

    U_a = -1.0;
    U_b =  2.0 * tan(Psi - 0.5 * M_PI) - 1.0;
}

///----------------------------------------------------------------------------

void  triag_st_make_simplex(double chi[3], const double eta, const double xi) noexcept {

    const double Xi_p_dsqrt3 =  xi / sqrt(3.0);

    chi[0] = 0.5 * (1.0 -  eta - Xi_p_dsqrt3 );
    chi[1] = 0.5 * (1.0 +  eta - Xi_p_dsqrt3 );
    chi[2] = Xi_p_dsqrt3;
}

/////////////////////////////////////////////////////////////////////////////////

template <>
void triag_ea_theta_limits<0>(double & theta_A, double & theta_B) noexcept {
    theta_A = 0.0;
    theta_B = M_PI / 3.0;
}

template <>
void triag_ea_theta_limits<1>(double & theta_A, double & theta_B) noexcept {
    theta_A = M_PI / 3.0;
    theta_B = M_PI / 2.0;
}

template <>
void triag_ea_theta_limits<2>(double & theta_A, double & theta_B) noexcept {
    theta_A = M_PI / 3.0;
    theta_B = M_PI / 2.0;
}

template <>
void triag_ea_theta_limits<3>(double & theta_A, double & theta_B) noexcept {
    theta_A = M_PI / 2.0;
    theta_B = M_PI;
}

template <>
void triag_ea_theta_limits<4>(double & theta_A, double & theta_B) noexcept {
    theta_A = 0.0;
    theta_B = M_PI / 2.0;
}

template <>
void triag_ea_theta_limits<5>(double & theta_A, double & theta_B) noexcept {

    theta_A = M_PI / 2.0;
    theta_B = M_PI;
}

///----------------------------------------------------------------------------

template <>
void quad_ea_theta_limits<0>(double & theta_A, double & theta_B) noexcept {
    theta_A = 0.0;
    theta_B = M_PI / 4.0;
}

template <>
void quad_ea_theta_limits<1>(double & theta_A, double & theta_B) noexcept {
    theta_A = 0.0;
    theta_B = M_PI / 4.0;
}

template <>
void quad_ea_theta_limits<2>(double & theta_A, double & theta_B) noexcept {
    theta_A = M_PI / 4.0;
    theta_B = M_PI / 2.0;
}

template <>
void quad_ea_theta_limits<3>(double & theta_A, double & theta_B) noexcept {
    theta_A = M_PI / 4.0;
    theta_B = M_PI / 2.0;
}

template <>
void quad_ea_theta_limits<4>(double & theta_A, double & theta_B) noexcept {
    theta_A = M_PI / 2.0;
    theta_B = 3.0 * M_PI / 4.0;
}

template <>
void quad_ea_theta_limits<5>(double & theta_A, double & theta_B) noexcept {
    theta_A = M_PI / 2.0;
    theta_B = 3.0 * M_PI / 4.0;
}

template <>
void quad_ea_theta_limits<6>(double & theta_A, double & theta_B) noexcept {
    theta_A = 3.0 * M_PI / 4.0;
    theta_B = M_PI;
}

template <>
void quad_ea_theta_limits<7>(double & theta_A, double & theta_B) noexcept {
    theta_A = 3.0 * M_PI / 4.0;
    theta_B = M_PI;
}

///----------------------------------------------------------------------------

template <>
void triag_ea_psi_limits<0>(double & psi_a, double & ) noexcept {
    psi_a = 0.0;
    //psi_b = psi_b;
}

template <>
void triag_ea_psi_limits<1>(double & , double & ) noexcept {
//    Default values
//    psi_a = psi_a;
//    psi_b = psi_b;
}

template <>
void triag_ea_psi_limits<2>(double & psi_a, double & psi_b) noexcept {
    psi_b = psi_a;
    psi_a = 0.0;
}

template <>
void triag_ea_psi_limits<3>(double & psi_a, double & psi_b) noexcept {
    psi_b = psi_a;
    psi_a = 0.0;
}

template <>
void triag_ea_psi_limits<4>(double & psi_a, double & psi_b) noexcept {
    psi_a = psi_b;
    psi_b = M_PI / 2.0;
}

template <>
void triag_ea_psi_limits<5>(double & , double & psi_b) noexcept {
    //psi_a = psi_a
    psi_b = M_PI / 2.0;
}

void  triag_ea_eta_limits(double & ua, double & ub) noexcept {
    ua = 0.0;
    ub = 1.0;
}

///----------------------------------------------------------------------------

template <>
void quad_ea_psi_limits<0>(double & psi_a, double & psi_b) noexcept {
    psi_b = psi_a;
    psi_a = 0.0;
}

template <>
void quad_ea_psi_limits<1>(double & , double & psi_b) noexcept {
//    psi_a = psi_a;
    psi_b = M_PI / 2.0;
}

template <>
void quad_ea_psi_limits<2>(double & psi_a, double & ) noexcept {
    psi_a = 0.0;
    //psi_b = psi_b;
}

template <>
void quad_ea_psi_limits<3>(double & psi_a, double & psi_b) noexcept {
    psi_a = psi_b;
    psi_b = M_PI / 2.0;
}

template <>
void quad_ea_psi_limits<4>(double & psi_a, double & ) noexcept {
    psi_a = 0.0;
    //psi_b = psi_b;
}

template <>
void quad_ea_psi_limits<5>(double & psi_a, double & psi_b) noexcept {
    psi_a = psi_b;
    psi_b = M_PI / 2.0;
}

template <>
void quad_ea_psi_limits<6>(double & psi_a, double & psi_b) noexcept {
    psi_b = -psi_a;
    psi_a = 0.0;
}

template <>
void quad_ea_psi_limits<7>(double & psi_a, double & psi_b) noexcept {
    psi_a *= -1.0;
    psi_b  = M_PI / 2.0;
}


void  quadr_ea_eta_limits(double & ua, double & ub) noexcept {
    ua = -1.0;
    ub =  1.0;
}

/////////////////////////////////////////////////////////////////////////////////

void triag_va_theta_p_limits(double & theta_p_A, double & theta_p_B) noexcept {

    theta_p_A = 0.0;
    theta_p_B = M_PI / 3.0;
}

void triag_va_theta_q_limits(double & theta_q_A, double & theta_q_B) noexcept {

    theta_q_A = 0.0;
    theta_q_B = M_PI / 3.0;
}

template <>
void quad_va_theta_p_limits<0>(double & theta_A, double & theta_B) noexcept {

    theta_A = 0.0;
    theta_B = M_PI / 4.0;
}

template <>
void quad_va_theta_p_limits<1>(double & theta_A, double & theta_B) noexcept {

    theta_A = 0.0;
    theta_B = M_PI / 4.0;
}

template <>
void quad_va_theta_p_limits<2>(double & theta_A, double & theta_B) noexcept {

    theta_A = M_PI / 4.0;
    theta_B = M_PI / 2.0;
}

template <>
void quad_va_theta_p_limits<3>(double & theta_A, double & theta_B) noexcept {

    theta_A = M_PI / 4.0;
    theta_B = M_PI / 2.0;
}

template <>
void quad_va_theta_q_limits<0>(double & theta_A, double & theta_B) noexcept {

    theta_A = 0.0;
    theta_B = M_PI / 4.0;
}

template <>
void quad_va_theta_q_limits<1>(double & theta_A, double & theta_B) noexcept {

    theta_A = M_PI / 4.0;
    theta_B = M_PI / 2.0;
}

template <>
void quad_va_theta_q_limits<2>(double & theta_A, double & theta_B) noexcept {

    theta_A = 0.0;
    theta_B = M_PI / 4.0;
}

template <>
void quad_va_theta_q_limits<3>(double & theta_A, double & theta_B) noexcept {

    theta_A = M_PI / 4.0;
    theta_B = M_PI / 2.0;
}

static double triangular_Lx_p_or_q(const double cosx_inp, const double sinx_inp) noexcept {

    const double sqrt3 = sqrt(3.0);
    return  2.0 * sqrt3 / ( sinx_inp + sqrt3 * cosx_inp);
}

double triag_va_Lp1(const double cosx_inp, const double sinx_inp) noexcept {
    return triangular_Lx_p_or_q(cosx_inp, sinx_inp);
}

double triag_va_Lq2(const double cosx_inp, const double sinx_inp) noexcept {
    return triangular_Lx_p_or_q(cosx_inp, sinx_inp);
}

template <>
double  quad_va_Lp1<0>(const double cosx_inp, const double ) noexcept {
    return  2.0 / cosx_inp;
}

template <>
double  quad_va_Lp1<1>(const double cosx_inp, const double ) noexcept {
    return  2.0 / cosx_inp;
}

template <>
double  quad_va_Lp1<2>(const double , const double sinx_inp) noexcept {
    return  2.0 / sinx_inp;
}

template <>
double quad_va_Lp1<3>(const double , const double sinx_inp) noexcept {
    return  2.0 / sinx_inp;
}

template <>
double quad_va_Lq2<0>(const double cosx_inp, const double ) noexcept {
    return 2.0 / cosx_inp;
}

template <>
double quad_va_Lq2<1>(const double , const double sinx_inp) noexcept {
    return 2.0 / sinx_inp;
}

template <>
double quad_va_Lq2<2>(const double cosx_inp, const double ) noexcept {
    return 2.0 / cosx_inp;
}

template <>
double quad_va_Lq2<3>(const double , const double sinx_inp) noexcept {
    return 2.0 / sinx_inp;
}

///////////////////////////////////////////////////////////////////////////////

} // End of the namespace Directfn

// End of the file

