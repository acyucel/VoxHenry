#ifndef _VERTEXADJACENT_ANYLATERAL_H_
#define _VERTEXADJACENT_ANYLATERAL_H_

#include "directfn_interface.h"
#include "directfn_common.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class DirectfnAlgorithm_VA : public virtual DirectfnInterface<ParticularKernel> {
public:

    DirectfnAlgorithm_VA();
    ~DirectfnAlgorithm_VA();

    DirectfnAlgorithm_VA(const DirectfnAlgorithm_VA & ) = delete;
    DirectfnAlgorithm_VA(DirectfnAlgorithm_VA && ) = delete;
    DirectfnAlgorithm_VA & operator = (const DirectfnAlgorithm_VA & ) = delete;
    DirectfnAlgorithm_VA & operator = (DirectfnAlgorithm_VA && ) = delete;

protected:
    /*! Pointer to the current theta_p_lim_1 */
    pFUN_lims_Ni    pf_theta_p_lim_1_crnt_;
    /*! Pointer to the current theta_q_lim_2 */
    pFUN_lims_Ni    pf_theta_q_lim_2_crnt_;
    /*! Pointer to the current Lp1 */
    pFUN_LpLq       pf_get_Lp_1_crnt_;
    /*! Pointer to the current Lq1*/
    pFUN_LpLq       pf_get_Lq_2_crnt_;
    /*! Pointer to the current lam max lim */
    pMFun_lam4_max<ParticularKernel>  pf_crnt_lam_maxlim_;

    /*! Memory for internal subintegrals */
    unique_ptr<dcomplex []>  Isub_;

private:
    /*! Prehashing constants. */
    double sin_Theta_p_1_;
    double cos_Theta_p_1_;
    double sin_Theta_q_2_;
    double cos_Theta_q_2_;
    double psi_a_;
    double psi_b_;
    double sin_Psi_3_;
    double cos_Psi_3_;
    double Lp_1_;
    double Lq_2_;

    virtual void  do_I_surface_surface_();
    void          gather_Iss_() noexcept;

    void     update_theta_p1_trigonometry_(const double Theta_p) noexcept;
    void     update_theta_q2_trigonometry_(const double Theta_q) noexcept;
    void     define_psi_a_psi_b_lamlim3_case_I_1_ (const double) noexcept;
    void     define_psi_a_psi_b_lamlim3_case_II_2_(const double) noexcept;
    void     calc_I_psi_N3_() noexcept;
    void     update_psi_3_trigonometry_(const double Psi_3) noexcept;
    void     calc_I_rho_4_() noexcept;
    double   lam_lim_b_1_() const noexcept;
    double   lam_lim_b_2_() const noexcept;
    void     calc_pq_simplex_(double xi_p_out[3], double xi_q_out[3], const double t_Lam) noexcept;

    virtual size_t   sub_ranges_numb_m_()  const noexcept = 0;
    virtual void     update_theta_p_lims_N1_fptr_(const size_t m) noexcept = 0;
    virtual void     update_theta_q_lims_N2_fptr_(const size_t m) noexcept = 0;
    virtual void     update_Lp1_fptr_(const size_t m) noexcept = 0;
    virtual void     update_Lq2_fptr_(const size_t m) noexcept = 0;
    virtual double   unity4_zero3_() const noexcept = 0;
    virtual void     va_make_simplex(double uvxi_p_out[3], const double u_eta_p,
                                             const double v_xi_p) const noexcept = 0;
};

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class Triangular_VA final :  public DirectfnAlgorithm_VA<ParticularKernel> {
public:
    Triangular_VA();
    ~Triangular_VA();

    Triangular_VA(const Triangular_VA & ) = delete;
    Triangular_VA(Triangular_VA && ) = delete;
    Triangular_VA & operator = (const Triangular_VA & ) = delete;
    Triangular_VA & operator = (Triangular_VA && ) = delete;

    virtual string name() const noexcept;

private:
    static constexpr size_t sub_figures_1_() noexcept {return 1;}

    void  allocate_I_vars_();
    void  initialize_limit_fptrs_() noexcept;

    virtual size_t   sub_ranges_numb_m_() const noexcept;
    virtual void     update_theta_p_lims_N1_fptr_(const size_t m) noexcept;
    virtual void     update_theta_q_lims_N2_fptr_(const size_t m) noexcept;
    virtual void     update_Lp1_fptr_(const size_t m) noexcept;
    virtual void     update_Lq2_fptr_(const size_t m) noexcept;
    virtual double   unity4_zero3_() const noexcept;
    virtual void     va_make_simplex(double uvxi_p_out[3], const double u_eta_p,
                                                           const double v_xi_p) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class Quadrilateral_VA : public DirectfnAlgorithm_VA<ParticularKernel> {
public:
    Quadrilateral_VA();
    ~Quadrilateral_VA();

    Quadrilateral_VA(const Quadrilateral_VA & ) = delete;
    Quadrilateral_VA(Quadrilateral_VA && ) = delete;
    Quadrilateral_VA & operator = (const Quadrilateral_VA & ) = delete;
    Quadrilateral_VA & operator = (Quadrilateral_VA && ) = delete;

    virtual string name() const noexcept;

protected:
    void  allocate_I_vars_();
    void  initialize_limit_fptrs_() noexcept;

private:
    /*! Forward declaration. Clang (at present) does not support constexpr function call. */
    static constexpr const size_t sub_figures_4_cval_ = 4;
    static constexpr size_t sub_figures_4_() noexcept {return sub_figures_4_cval_;}

    /*! Array of length 6 of pointers to: psi_limits(m, &PSI_A, &PSI_B); m = 1,2,... 6  */
    pFUN_lims_Ni  pF_theta_p_lim1_arr_[sub_figures_4_cval_];

    /*! Array of length 6 of pointers to: lam_limits(m, PSI, ETA, &LAM_A, &LAM_B); m = 1,2,... 6 */
    pFUN_lims_Ni  pF_theta_q_lim2_arr_[sub_figures_4_cval_];

    pFUN_LpLq  pf_Lp_1_arr_[sub_figures_4_cval_];
    pFUN_LpLq  pf_Lq_2_arr_[sub_figures_4_cval_];

    virtual size_t   sub_ranges_numb_m_()  const noexcept;
    virtual void     update_theta_p_lims_N1_fptr_(const size_t m) noexcept;
    virtual void     update_theta_q_lims_N2_fptr_(const size_t m) noexcept;
    virtual void     update_Lp1_fptr_(const size_t m) noexcept;
    virtual void     update_Lq2_fptr_(const size_t m) noexcept;
    virtual double   unity4_zero3_() const noexcept;
    virtual void     va_make_simplex(double uvxi_p_out[3], const double u_eta_p,
                                                           const double v_xi_p) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  // _VERTEXADJACENT_ANYLATERAL_H_

// End of the file

