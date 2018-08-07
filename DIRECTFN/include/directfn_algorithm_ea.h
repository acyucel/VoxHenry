#ifndef _EDGEADJACENT_ANYLATERAL_H_
#define _EDGEADJACENT_ANYLATERAL_H_

#include "directfn_interface.h"
#include "directfn_common.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class DirectfnAlgorithm_EA : public virtual DirectfnInterface<ParticularKernel> {
public:
    DirectfnAlgorithm_EA();
    virtual ~DirectfnAlgorithm_EA();

    DirectfnAlgorithm_EA(const DirectfnAlgorithm_EA & ) = delete;
    DirectfnAlgorithm_EA(DirectfnAlgorithm_EA && ) = delete;
    DirectfnAlgorithm_EA & operator = (const DirectfnAlgorithm_EA & ) = delete;
    DirectfnAlgorithm_EA & operator = (DirectfnAlgorithm_EA && ) = delete;

protected:
    /*! Function pointer to the first limits for N1 */
    pFUN_ea_theta_lim_N1       pf_thetalim_crnt_;
    pFUN_ea_psi_lim_N1         pf_psilim_crnt_;
    /*! Current pointer to the m-cases for N3 */
    pMFUN_ea_alg_nfun_IorII<ParticularKernel>    fptr_nfunc_crnt_;

    double  cos_theta_1_;
    double  sin_theta_1_;
    double  tan_theta_1_;
    double  tan_05pi_mns_theta_1_;
    double  cos_psi_2_;
    double  sin_psi_2_;
    double  tan_psi_2_;

    /*! Memory for internal subintegrals */
    unique_ptr<dcomplex []>  Isub_;

    /*! Declaration for reference to member function pointers
     *  to parent class from inherited classes: */
    friend class Triangular_EA<ParticularKernel>;
//    friend class Quadrilateral_EA<ParticularKernel>;

    /*! n_function with one call of a_function. */
    void  single_n_func_I_(const double mlt_factor) noexcept;
    /*! n_function with two calls of a_function. */
    void  double_n_func_II_(const double mlt_factor) noexcept;

private:
    /*! The implementation of integration algorithm for
     *  triangular and quadrilateral EA term.  */
    virtual void   do_I_surface_surface_();

    void  gather_Iss_() noexcept;
    void  update_theta_trigonometry_(const double Theta) noexcept;
    void  update_psi_trigonometry_(const double Psi) noexcept;

    void  calc_n_func_() noexcept;
    void  calc_afuncs_(const double Eta_3) noexcept;

    virtual size_t  sub_ranges_numb_m_() const noexcept = 0;
    virtual void    update_theta_lims_N1_fptr_(const size_t m) noexcept = 0;
    virtual void    update_psi_lims_N2_fptr_(const size_t m) noexcept = 0;
    virtual void    update_eta_lims_N3_fptr_(const size_t m) noexcept = 0;
    virtual void    update_lam_lims_N4_fprt_(const size_t m) noexcept = 0;

    virtual double  init_psi_a_by_theta_() noexcept = 0;
    virtual double  init_psi_b_by_theta_() noexcept = 0;
    virtual void    calc_I_lam_N3_nfuncs_wrapper_(const double mlt_factor) noexcept = 0;

    virtual void    setup_limits34_pntr_to_single_() noexcept = 0;
    virtual void    setup_limits34_pntr_to_double_first_() noexcept = 0;
    virtual void    setup_limits34_pntr_to_double_secnd_() noexcept = 0;

    virtual void    get_ueta_limits3_nfun_(double & , double & ) noexcept = 0;
    virtual double  get_lam_limit_(const double Eta_3) noexcept = 0;

    virtual void    ea_calc_pq_simplex_(double uvxi_p[3], double uvxi_q[3],
                                        const double Eta_3, const double Lam_4) noexcept = 0;
};

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class Triangular_EA : public virtual DirectfnAlgorithm_EA<ParticularKernel> {
public:

    Triangular_EA();
    virtual ~Triangular_EA();

    Triangular_EA(const Triangular_EA & ) = delete;
    Triangular_EA(Triangular_EA && ) = delete;
    Triangular_EA & operator = (const Triangular_EA & ) = delete;
    Triangular_EA & operator = (Triangular_EA && ) = delete;

    virtual string name() const noexcept;

private:
    /*! Forward declaration */
    static constexpr size_t sub_ranges_6_()  noexcept {return 6;}

    /*! Array of pointers to theta_1 limits for N1 integration. */
    pFUN_ea_theta_lim_N1      pf_thetalim_arr_[sub_ranges_6_()];
    pFUN_ea_psi_lim_N1        pf_psilim_arr_[sub_ranges_6_()];
    /*! Pointer array to the m-cases for N3 */
    pMFUN_ea_alg_nfun_IorII<ParticularKernel>   fptr_nfunc_arr_[sub_ranges_6_()];

    /*! This factor is used for simplex in the N4 loop.
     *  It is reset twice: before each call of n_functions */
    double eta_pq_sgn_;

    /*! Eta limits depending on single/double_IorII choice */
    pmFUN_ea_eta3_lim_N3_TR<ParticularKernel>  pf_eta3_lim_crnt_;
    pmFUN_ea_eta3_lim_N3_TR<ParticularKernel>  pf_eta3_lim_single_crnt_;
    pmFUN_ea_eta3_lim_N3_TR<ParticularKernel>  pf_eta3_lim_double_frst_crnt_;
    pmFUN_ea_eta3_lim_N3_TR<ParticularKernel>  pf_eta3_lim_double_scnd_crnt_;

    pmFUN_ea_eta3_lim_N3_TR<ParticularKernel>  pf_eta3_lim_single_arr_[sub_ranges_6_()];
    pmFUN_ea_eta3_lim_N3_TR<ParticularKernel>  pf_eta3_lim_double_frst_arr_[sub_ranges_6_()];
    pmFUN_ea_eta3_lim_N3_TR<ParticularKernel>  pf_eta3_lim_double_scnd_arr_[sub_ranges_6_()];

    /*! Lambda limits depending on single/double_IorII choice */
    pmFUN_ea_lam4_lim_N4_TR<ParticularKernel>  pf_lam4_lim_crnt_;
    pmFUN_ea_lam4_lim_N4_TR<ParticularKernel>  pf_lam4_lim_single_crnt_;
    pmFUN_ea_lam4_lim_N4_TR<ParticularKernel>  pf_lam4_lim_double_frst_crnt_;
    pmFUN_ea_lam4_lim_N4_TR<ParticularKernel>  pf_lam4_lim_double_scnd_crnt_;

    pmFUN_ea_lam4_lim_N4_TR<ParticularKernel>  pf_lam4_lim_single_arr_[sub_ranges_6_()];
    pmFUN_ea_lam4_lim_N4_TR<ParticularKernel>  pf_lam4_lim_double_frst_arr_[sub_ranges_6_()];
    pmFUN_ea_lam4_lim_N4_TR<ParticularKernel>  pf_lam4_lim_double_scnd_arr_[sub_ranges_6_()];


    void    allocate_I_vars_() noexcept;
    void    initialize_limit_fptrs_() noexcept;

    void    eta3_lims_single_(double & Ua, double & Ub) noexcept;
    void    eta3_lims_double_frst_I_12(double & Ua, double & Ub) noexcept;
    void    eta3_lims_double_frst_I_3(double & Ua, double & Ub) noexcept;
    void    eta3_lims_double_scnd_II_12_(double & Ua, double & Ub) noexcept;
    void    eta3_lims_double_scnd_II_3_(double & Ua, double & Ub) noexcept;

    double  gamma_sin_cos_() const noexcept;
    double  gamma_tan_() const noexcept;

    double  lambda_L1_(const double ) const noexcept;
    double  lambda_L2_(const double ) const noexcept;
    double  lambda_L2_abs_(const double ) const noexcept;
    double  lambda_L3_(const double ) const noexcept;

    virtual size_t   sub_ranges_numb_m_() const noexcept;
    virtual void     update_theta_lims_N1_fptr_(const size_t m) noexcept;
    virtual void     update_psi_lims_N2_fptr_(const size_t m) noexcept;
    virtual void     update_eta_lims_N3_fptr_(const size_t m) noexcept;
    virtual void     update_lam_lims_N4_fprt_(const size_t m) noexcept;

    virtual double   init_psi_a_by_theta_() noexcept;
    virtual double   init_psi_b_by_theta_() noexcept;
    virtual void     calc_I_lam_N3_nfuncs_wrapper_(const double mlt_factor) noexcept;

    virtual void     setup_limits34_pntr_to_single_() noexcept;
    virtual void     setup_limits34_pntr_to_double_first_() noexcept;
    virtual void     setup_limits34_pntr_to_double_secnd_() noexcept;
    virtual void     get_ueta_limits3_nfun_(double & , double & ) noexcept;
    virtual double   get_lam_limit_(const double eta_3) noexcept;
    virtual void     ea_calc_pq_simplex_(double uvxi_p[3], double uvxi_q[3],
                                         const double Eta_3, const double Lam_4) noexcept;
};

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class Quadrilateral_EA : public DirectfnAlgorithm_EA<ParticularKernel> {
public:
    Quadrilateral_EA();
    virtual ~Quadrilateral_EA();

    Quadrilateral_EA(const Quadrilateral_EA & ) = delete;
    Quadrilateral_EA(Quadrilateral_EA && ) = delete;
    Quadrilateral_EA & operator = (const Quadrilateral_EA & ) = delete;
    Quadrilateral_EA & operator = (Quadrilateral_EA && ) = delete;

    virtual string name() const noexcept;

protected:
    void  allocate_I_vars_() noexcept;
    void  initialize_limit_fptrs_() noexcept;

private:
    static constexpr const size_t sub_ranges_8_cval_  = 8;
    /*! Forward declaration. */
    static constexpr size_t sub_ranges_8_()  noexcept {return sub_ranges_8_cval_;}

    /*! Array of pointers to theta_1 limits for N1 integration. */
    pFUN_ea_theta_lim_N1      pf_thetalim_arr_[sub_ranges_8_cval_];
    pFUN_ea_psi_lim_N1        pf_psilim_arr_[sub_ranges_8_cval_];
    /*! Pointer array to the m-cases for N3 */
    pMFUN_ea_alg_nfun_IorII<ParticularKernel>   fptr_nfunc_arr_[sub_ranges_8_cval_];

    /*! Eta limits depending on single/double_IorII choice */
    pmFUN_ea_eta3_lim_N3_QD<ParticularKernel>  pf_eta3_lim_crnt_;
    pmFUN_ea_eta3_lim_N3_QD<ParticularKernel>  pf_eta3_lim_single_crnt_;
    pmFUN_ea_eta3_lim_N3_QD<ParticularKernel>  pf_eta3_lim_double_frst_crnt_;
    pmFUN_ea_eta3_lim_N3_QD<ParticularKernel>  pf_eta3_lim_double_scnd_crnt_;

    pmFUN_ea_eta3_lim_N3_QD<ParticularKernel>  pf_eta3_lim_single_arr_[sub_ranges_8_cval_];
    pmFUN_ea_eta3_lim_N3_QD<ParticularKernel>  pf_eta3_lim_double_frst_arr_[sub_ranges_8_cval_];
    pmFUN_ea_eta3_lim_N3_QD<ParticularKernel>  pf_eta3_lim_double_scnd_arr_[sub_ranges_8_cval_];

    /*! Lambda limits depending on single/double_IorII choice */
    pmFUN_ea_lam4_lim_N4_QD<ParticularKernel>  pf_lam4_lim_crnt_;
    pmFUN_ea_lam4_lim_N4_QD<ParticularKernel>  pf_lam4_lim_single_crnt_;
    pmFUN_ea_lam4_lim_N4_QD<ParticularKernel>  pf_lam4_lim_double_frst_crnt_;
    pmFUN_ea_lam4_lim_N4_QD<ParticularKernel>  pf_lam4_lim_double_scnd_crnt_;

    pmFUN_ea_lam4_lim_N4_QD<ParticularKernel>  pf_lam4_lim_single_arr_[sub_ranges_8_cval_];
    pmFUN_ea_lam4_lim_N4_QD<ParticularKernel>  pf_lam4_lim_double_frst_arr_[sub_ranges_8_cval_];
    pmFUN_ea_lam4_lim_N4_QD<ParticularKernel>  pf_lam4_lim_double_scnd_arr_[sub_ranges_8_cval_];


    void  eta3_lims_single_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_m1Fm_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_m1Tm_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_m1Fp_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_m1Tp_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_Fmp1_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_Tmp1_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_Fpp1_(double & Ua, double & Ub) noexcept;
    void  eta3_lims_double_Tpp1_(double & Ua, double & Ub) noexcept;

    double F_plus_1_() const noexcept;
    double F_mnus_1_() const noexcept;
    double T_plus_1_() const noexcept;
    double T_mnus_1_() const noexcept;

    double  Lam_xp1_cc_(const double ) const noexcept;
    double  Lam_xm1_cc_(const double ) const noexcept;
    double  Lam_0p2_sc_(const double ) const noexcept;
    double  Lam_0p2_0s_(const double ) const noexcept;

    virtual size_t   sub_ranges_numb_m_() const noexcept;
    virtual void     update_theta_lims_N1_fptr_(const size_t m) noexcept;
    virtual void     update_psi_lims_N2_fptr_(const size_t m) noexcept;
    virtual void     update_eta_lims_N3_fptr_(const size_t m) noexcept;
    virtual void     update_lam_lims_N4_fprt_(const size_t m) noexcept;

    virtual double   init_psi_a_by_theta_() noexcept;
    virtual double   init_psi_b_by_theta_() noexcept;
    virtual void     calc_I_lam_N3_nfuncs_wrapper_(const double mlt_factor) noexcept;

    virtual void    setup_limits34_pntr_to_single_() noexcept;
    virtual void    setup_limits34_pntr_to_double_first_() noexcept;
    virtual void    setup_limits34_pntr_to_double_secnd_() noexcept;
    virtual void    get_ueta_limits3_nfun_(double & , double & ) noexcept;
    virtual double  get_lam_limit_(const double eta_3) noexcept;
    virtual void    ea_calc_pq_simplex_(double uvxi_p[3], double uvxi_q[3],
                                        const double Eta_3, const double Lam_4) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif   // _EDGEADJACENT_ANYLATERAL_H_

// End of the file

