#ifndef _SELFTERM_ANYLATERAL_H_
#define _SELFTERM_ANYLATERAL_H_

#include "directfn_interface.h"
#include "directfn_common.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class DirectfnAlgorithm_ST : public virtual DirectfnInterface<ParticularKernel> {
public:
    DirectfnAlgorithm_ST();
    virtual ~DirectfnAlgorithm_ST();

    DirectfnAlgorithm_ST(const DirectfnAlgorithm_ST & ) = delete;
    DirectfnAlgorithm_ST(DirectfnAlgorithm_ST && ) = delete;
    DirectfnAlgorithm_ST & operator = (const DirectfnAlgorithm_ST & ) = delete;
    DirectfnAlgorithm_ST & operator = (DirectfnAlgorithm_ST && ) = delete;

protected:

    /*! Hashed variable(s) precomputed in the upper loop(s) */
    double sin_PSI_;
    double cos_PSI_;
    double t_Eta_2_;

    /*! Memory for internal subintegrals */
    unique_ptr<dcomplex []>  Isub_;
    unique_ptr<dcomplex []>  Isd_k_;

    /*! A pointer to a particular function from the  pf_subtrg_arr_ array.
     *  Is called from a_functions( PSI,ETA,LAMBDA,kk);  kk = 1, 2, 3, (4) */
    pFUN_tq_st_SUBTR  pf_subtrg_crnt_;

    /*! Pointer to the current psi_lim */
    pFUN_st_psilim_N1  pf_psilim_crnt_;

    /*! Pointer to the current psi_lim */
    pFUN_st_etalim_N2  pf_etalim_crnt_;

    virtual void   precompute_psi_trigonomety_(const double t_PSI) noexcept;

private:
    /*! The implementation of integration algorithm for
     *  triangular and quadrilateral ST term.  */
    virtual void   do_I_surface_surface_();

    /*! Auxiliary routines for internal intermediate integrals save */
    void    init_Isub_() noexcept;
    void    gather_Isub_(const size_t kk) noexcept;
    void    gather_Iss_() noexcept;

    virtual size_t   sub_figures_numb_k_() const noexcept = 0;
    virtual size_t   sub_ranges_numb_m_()  const noexcept = 0;

    virtual void     update_subtri_fptrs_(const size_t kk) noexcept = 0;
    virtual void     update_psilims_fptrs_N1_(const size_t m) noexcept = 0;
    virtual void     update_etalims_fptrs_N2_(const size_t m) noexcept = 0;
    virtual void     update_lamlims_fptrs_N3_(const size_t m) noexcept = 0;
    virtual double   tri_zero_qua_one_() const noexcept = 0;
    virtual void     st_make_simplex_(double chi[3], const double eta, const double xi) noexcept = 0;
    virtual void     setup_lamlim_(const double Eta_2, double & Lam_b) const noexcept = 0;

    /*! N3 loop is done inside. Vectorized summation over kernel index-i
     *  is applied inside the n_functions for lam_3. */
    void calc_n_functions_() noexcept;

    /*! Peform the forth integration for N4 by the Rho.
     *  kk and m arguments have been setup in main loop before enter in this call.
     *  Vectorized summation over kernel index-i
     *  is applied inside the a_functions for rho_4*/
    void calc_a_functions_(const double t_Lam) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class Triangular_ST final : public DirectfnAlgorithm_ST<ParticularKernel> {
public:
    Triangular_ST();
    ~Triangular_ST();

    Triangular_ST(const Triangular_ST & ) = delete;
    Triangular_ST(Triangular_ST && ) = delete;
    Triangular_ST & operator = (const Triangular_ST & ) = delete;
    Triangular_ST & operator = (Triangular_ST && ) = delete;

    virtual string name() const noexcept;
private:
    /*! Forward declaration */
    static constexpr size_t sub_ranges_8_()  noexcept {return 8;}
    static constexpr size_t sub_figures_3_() noexcept {return 3;}

    /*! Array of pointers to templated functions. Can be used in loops as
     *  a reference to a particular template,  when template's index
     *  is not constexpr (needed for templates specialization)  */
    pFUN_tq_st_SUBTR  pf_subtrg_arr_[sub_figures_3_()];

    /*! Array of length 8 of pointers to: psi_limits(m, &PSI_A, &PSI_B); m = 1,2,... 8  */
    pFUN_st_psilim_N1  pf_psilim_arr_[sub_ranges_8_()];

    /*! Array of length 8 of pointers to: eta_limits(m, &PSI_A, &PSI_B); m = 1,2,... 8  */
    pFUN_st_etalim_N2  pf_etalim_arr_[sub_ranges_8_()];

    /*! Array of length 8 of pointers to: lam_limits(m, PSI, ETA, &LAM_A, &LAM_B); m = 1,2,... 8  */
    pMFUN_tr_st_lamlim_N3<ParticularKernel>  pmf_lamlim_arr_[sub_ranges_8_()];

    /*! Pointer to the current eta_lim */
    pMFUN_tr_st_lamlim_N3<ParticularKernel>  pmf_lamlim_crnt_;

    void  allocate_I_vars_();
    void  initialize_limit_fptrs_();

    virtual size_t  sub_figures_numb_k_() const noexcept;
    virtual size_t  sub_ranges_numb_m_()  const noexcept;

    virtual void    update_subtri_fptrs_(const size_t kk) noexcept;
    virtual void    update_psilims_fptrs_N1_(const size_t m) noexcept;
    virtual void    update_etalims_fptrs_N2_(const size_t m) noexcept;
    virtual void    update_lamlims_fptrs_N3_(const size_t m) noexcept;
    virtual double  tri_zero_qua_one_() const noexcept;
    virtual void    st_make_simplex_(double chi[3], const double eta, const double xi) noexcept;

    virtual void    setup_lamlim_(const double Eta_2, double & Lam_b) const noexcept;

    void  triag_st_lam_limits_1_(const double Eta, double & Lam_b) const noexcept;
    void  triag_st_lam_limits_2_(const double Eta, double & Lam_b) const noexcept;
    void  triag_st_lam_limits_3_(const double Eta, double & Lam_b) const noexcept;
    void  triag_st_lam_limits_4_(const double Eta, double & Lam_b) const noexcept;
    void  triag_st_lam_limits_5_(const double Eta, double & Lam_b) const noexcept;
    void  triag_st_lam_limits_6_(const double Eta, double & Lam_b) const noexcept;
    void  triag_st_lam_limits_7_(const double Eta, double & Lam_b) const noexcept;
    void  triag_st_lam_limits_8_(const double Eta, double & Lam_b) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

/*! \class Quadrilateral_ST   */

template <typename ParticularKernel>
class Quadrilateral_ST : public DirectfnAlgorithm_ST<ParticularKernel> {
public:

    Quadrilateral_ST();
    virtual ~Quadrilateral_ST();

    Quadrilateral_ST(const Quadrilateral_ST & ) = delete;
    Quadrilateral_ST(Quadrilateral_ST && ) = delete;
    Quadrilateral_ST & operator = (const Quadrilateral_ST & ) = delete;
    Quadrilateral_ST & operator = (Quadrilateral_ST && ) = delete;

    /*! Particular name of the class. */
    virtual string name() const noexcept;

protected:
    void  allocate_I_vars_();
    void  allocate_I_vars_(const size_t ker_sz) noexcept;
    void  initialize_limit_fptrs_();

private:
    static constexpr const size_t sub_ranges_6_cval_  = 6;
    static constexpr const size_t sub_figures_4_cval_ = 4;

    /*! Forward declaration. */
    static constexpr size_t sub_ranges_6_()  noexcept {return sub_ranges_6_cval_;}
    static constexpr size_t sub_figures_4_() noexcept {return sub_figures_4_cval_;}

    /*! Precompured value in the N1 loop. */
    double cos_pi_ms_PSI_;

    /*! Array of pointers to templated functions. Can be used in loops as
     *  a reference to a particular template,  when template's index
     *  is not constexpr (needed for templates specialization)  */
    pFUN_tq_st_SUBTR   pf_subtrg_arr_[sub_figures_4_cval_]; //[sub_figures_4_()]; <-- clang bug, change when fixed

    /*! Array of length 6 of pointers to: psi_limits(m, &PSI_A, &PSI_B); m = 1,2,... 6  */
    pFUN_st_psilim_N1  pf_psilim_arr_[sub_ranges_6_cval_]; //[sub_ranges_6_()]; <-- clang bug, change when fixed

    /*! Array of length 6 of pointers to: lam_limits(m, PSI, ETA, &LAM_A, &LAM_B); m = 1,2,... 6 */
    pFUN_st_etalim_N2  pf_etalim_arr_[sub_ranges_6_cval_]; // [sub_ranges_6_()]; <-- clang bug, change when fixed

    /*! Array of length 6 of pointers to: lam_limits(m, PSI, ETA, &LAM_A, &LAM_B); m = 1,2,... 6  */
    pMFUN_qd_st_lamlim_N3<ParticularKernel>  pmf_lamlim_arr_[sub_ranges_6_cval_];  // [sub_ranges_6_()];

    /*! Pointer to the current eta_lim */
    pMFUN_qd_st_lamlim_N3<ParticularKernel>  pmf_lamlim_crnt_;

    virtual void   precompute_psi_trigonomety_(const double t_PSI) noexcept final;

    virtual size_t  sub_figures_numb_k_() const noexcept final;
    virtual size_t  sub_ranges_numb_m_()  const noexcept final;

    virtual void    update_subtri_fptrs_(const size_t kk) noexcept final;
    virtual void    update_psilims_fptrs_N1_(const size_t m) noexcept final;
    virtual void    update_etalims_fptrs_N2_(const size_t m) noexcept final;
    virtual void    update_lamlims_fptrs_N3_(const size_t m) noexcept final;
    virtual double  tri_zero_qua_one_() const noexcept final;
    virtual void    st_make_simplex_(double chi[3], const double eta, const double xi) noexcept final;
    virtual void    setup_lamlim_(const double Eta_2, double & Lam_b) const noexcept final;

    void  quad_st_lamlims_12_(const double U, double & LAM_b) const noexcept;
    void  quad_st_lamlims_34_(const double U, double & LAM_b) const noexcept;
    void  quad_st_lamlims_56_(const double U, double & LAM_b) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif   // _SELFTERM_ANYLATERAL_H_

// End of the file

