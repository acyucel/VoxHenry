#ifndef _DIRECTFN_COMMON_H_
#define _DIRECTFN_COMMON_H_

#include "directfn_defs.h"
//#include <cfloat>
//#include <memory>
//#include <limits>

//using std::numeric_limits<double>::epsilon();
namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

/*! Gauss quadrature points and weights  */
void gl_xw_1d (int n, double x[], double w[]);

/*! Computes the dot product x[0] * y[0] + x[1] * y[1] + x[2] * y[2].
 *  No complex conjugation for neither the first nor the second arguments.
 *
 * \param x - the first  vector of the length 3 of type double or complex.
 * \param y - the second vector of the length 3 of type double or complex.
 * \return scalar value of the double or complex<double> type
 */
template <typename T_first_input, typename T_secnd_input, typename T_result>
T_result  vector_dot(const T_first_input x[3], const T_secnd_input y[3]);
double    vector_dot(const double x[3], const double y[3]);

double    sqrt_vector_dot(const double x[3], const double y[3]);

/*! Computes the cross product of two length-3 vectors.
 *  \param[in]  x - the first  vector of the length 3 of type double or complex.
 *  \param[in]  y - the second vector of the length 3 of type double or complex.
 *  \param[out] z - result */
template <typename Real_T, typename Complex_T>
void vector_cross(const Complex_T x[3], const Real_T y[3], Complex_T z[3]);

/*!  Computes normal np to a given triangle defined by 3 vertexes
 *  \param[in]  r1 - first  vertex radius,
 *  \param[in]  r2 - second vertex radius,
 *  \param[in]  r3 - third  vertex radius,
 *  \param[out] nr - resulting normal.
 */
void  calc_normal(double nr[3], const double r1[3], const double r2[3], const double r3[3]) noexcept;

//class DirectfnInterface;
//void  print(const DirectfnInterface * , const dcomplex * ref_val, const dcomplex * eth_val) noexcept;

void calc_rotation_matrix_co(const double u[3], const double theta, double r[3]) noexcept;


//double max_element(double Error[]);
double max_element(const double * const Error, const size_t err_len);

#ifndef DBL_EPSILON
const double DBL_EPSILON = 3e-16;
#endif
double relative_error(dcomplex x, dcomplex x_ref);

///////////////////////////////////////////////////////////////////////////////
/// ST

using std::size_t;

/*! in, in, out, out */
using pFUN_tq_st_SUBTR = void (*) (const double , const double , double & , double & );

/*! Pointer type to the t/q_st_psi_limits function. */
using pFUN_st_psilim_N1 = void (*) (double & , double & ); // noexcept;

/*! Pointer type to the t/q_st_psi_limits function. */
using pFUN_st_etalim_N2 = void (*) (const double, double & , double &);


/*! Defines psi limits for N1 integration for triangular element ST */
template <size_t argument>
void  triag_st_psi_limits(double & psi_A, double & psi_B) noexcept;

/*! Defines psi limits for N1 integration for quadrilateral element ST */
template <size_t>
void  quad_st_psi_limits(double & psi_A, double & psi_B) noexcept;


///*! Defines eta limits for N2 integration for triangular element ST */
template <size_t argument>
void  triag_st_eta_limits(const double Psi, double & eta_A, double & eta_B) noexcept;

/*! Defines eta limits for N2 integration for quadrilateral element ST */
template <size_t argument>
void  quad_st_eta_limits(const double Psi, double & eta_A, double & eta_B) noexcept;


/*!  Forward declaration of templated classes. */
template <typename ParticularKernel> class Triangular_ST;
template <typename ParticularKernel> class Quadrilateral_ST;

///*! Defines lambda limites for N3 integration for triangles */
template <typename ParticularKernel>
using pMFUN_tr_st_lamlim_N3 = void (Triangular_ST<ParticularKernel>::*)(const double , double &) const;

/*! Defines lambda limites for N3 integration for quadrangles */
template <typename ParticularKernel>
using pMFUN_qd_st_lamlim_N3 = void (Quadrilateral_ST<ParticularKernel>::*)(const double , double &) const;


///*! Subtriangles for Triangular elements */
template <size_t argument>
void triag_st_subtriangles(const double ETA_sub, const double XI_sub, double & ETA, double & XI) noexcept;

/*! Subtriangles for Quadrilateral elements */
template <size_t argument>
void  quadr_st_subtriangles(const double U_sub, const double V_sub, double & U, double & V) noexcept;

/*!  Update xi_p. For q - call again. */
void  triag_st_make_simplex(double chi[3], const double eta, const double xi) noexcept;

///////////////////////////////////////////////////////////////////////////////

///----------------------------------------------------------------------------
/// EA

template <size_t  m_index>
void triag_ea_theta_limits(double & , double & ) noexcept;

template <size_t  m_index>
void quad_ea_theta_limits(double & , double & ) noexcept;

/*! Pointer type to the t/q_st_psi_limits function. */
using pFUN_ea_theta_lim_N1 = void (*) (double & , double &);

/*! The input param-1 and param-2 are initializers of psi_a and psi_b */
template <size_t m_index>
void  triag_ea_psi_limits(double & , double & ) noexcept;

/*! The input param-1 and param-2 are initializers of psi_a and psi_b */
template <size_t m_index>
void  quad_ea_psi_limits(double & , double & ) noexcept;

using pFUN_ea_psi_lim_N1 = void (*) (double & , double &);

/*! Forward declaration */
template <typename ParticularKernel> class DirectfnAlgorithm_EA;
template <typename ParticularKernel> class Triangular_EA;
template <typename ParticularKernel> class Quadrilateral_EA;

template <typename ParticularKernel>
using pMFUN_ea_alg_nfun_IorII = void (DirectfnAlgorithm_EA<ParticularKernel>::*)(const double mlt_fctr);

void  triag_ea_eta_limits(double & ua, double & ub) noexcept;
void  quadr_ea_eta_limits(double & ua, double & ub) noexcept;

template <typename ParticularKernel>
using pmFUN_ea_eta3_lim_N3_TR = void (Triangular_EA<ParticularKernel>::*)(double & , double &);

template <typename ParticularKernel>
using pmFUN_ea_lam4_lim_N4_TR = double (Triangular_EA<ParticularKernel>::*)(const double) const;

template <typename ParticularKernel>
using pmFUN_ea_eta3_lim_N3_QD = void (Quadrilateral_EA<ParticularKernel>::*) (double & , double &);

template <typename ParticularKernel>
using pmFUN_ea_lam4_lim_N4_QD = double (Quadrilateral_EA<ParticularKernel>::*) (const double) const;

/////////////////////////////////////////////////////////////////////////////////

///----------------------------------------------------------------------------
/// VA

///*! Defines limits for the N1 integration for triangular element VA */
void triag_va_theta_p_limits(double &, double & ) noexcept;

void triag_va_theta_q_limits(double &, double & ) noexcept;

/*! Defines limits for the N1 integration for quadrilateral element VA */
template <size_t m_index>
void quad_va_theta_p_limits(double &, double & ) noexcept;

template <size_t m_index>
void quad_va_theta_q_limits(double &, double & ) noexcept;

using pFUN_lims_Ni = void (*) (double & , double & );

double  triag_va_Lp1(const double cosx_inp, const double sinx_inp) noexcept;

double  triag_va_Lq2(const double cosx_inp, const double sinx_inp) noexcept;

template <size_t>
double  quad_va_Lp1(const double cosx_inp, const double sinx_inp) noexcept;

template <size_t>
double  quad_va_Lq2(const double cosx_inp, const double sinx_inp) noexcept;

using pFUN_LpLq = double (*) (const double , const double );

/*! Forward declaration */
template <typename ParticularKernel> class DirectfnAlgorithm_VA;

/*! Pointer to the member function */
template <typename ParticularKernel>
using pMFun_lam4_max = double (DirectfnAlgorithm_VA<ParticularKernel>::*) () const;

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif   // _DIRECTFN_COMMON_H_

// End of the file

