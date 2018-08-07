#include <iostream>
#include <iomanip>
#include "directfn_algorithm_ea.h"
#include "directfn_kernel_tri.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"
#include "directfn_kernel_quad_voxhenry.h"

using std::cout;
using std::endl;
using std::setprecision;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
DirectfnAlgorithm_EA<ParticularKernel>::DirectfnAlgorithm_EA():
DirectfnInterface<ParticularKernel>(),
pf_thetalim_crnt_(nullptr),
pf_psilim_crnt_(nullptr),
fptr_nfunc_crnt_(nullptr),
cos_theta_1_(0.0),
sin_theta_1_(0.0),
tan_05pi_mns_theta_1_(0.0),
cos_psi_2_(0.0),
sin_psi_2_(0.0),
tan_psi_2_(0.0),
Isub_(nullptr) {
}

//virtual
template <typename ParticularKernel>
DirectfnAlgorithm_EA<ParticularKernel>::~DirectfnAlgorithm_EA() {
}

//virtual
template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::do_I_surface_surface_() {

    for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {

        // Bind function pointers to get the limits of integration variables:
        update_theta_lims_N1_fptr_(m);
        update_psi_lims_N2_fptr_(m);
        update_eta_lims_N3_fptr_(m);  // Bind function pointers to single or double N3 afunc calls. Also bind limits.
        update_lam_lims_N4_fprt_(m);  // Update lambda_b pointers

        double theta_A = 0.0, theta_B = 0.0;
        pf_thetalim_crnt_(theta_A, theta_B);
        const double theta_a_plus_b = 0.5 * (theta_A + theta_B);
        const double theta_b_mnus_a = 0.5 * (theta_B - theta_A);
        const double * const cz1 = this->up_z1_.get();
        const double * const cw1 = this->up_w1_.get();

        this->up_kerSummator_->nullify_Ipsi_1();
        for (size_t n_theta_1 = 0; n_theta_1 < this->N1(); ++n_theta_1) {

            const double Theta = theta_b_mnus_a * cz1[n_theta_1] + theta_a_plus_b;

            update_theta_trigonometry_(Theta);
            double psi_a = init_psi_a_by_theta_();  // Theta dependent cos/sin have been already prehashed
            double psi_b = init_psi_b_by_theta_();  // Theta dependent cos/sin have been already prehashed
            // Depends on psi_a(theta), psi_b(theta) also, so it is called here, not above
            pf_psilim_crnt_(psi_a, psi_b); // redefini psi_a/b according to cases
            const double psi_b_mnus_a = 0.5 * (psi_b - psi_a);
            const double psi_a_plus_b = 0.5 * (psi_a + psi_b);
            const double * const cz2 = this->up_z2_.get();
            const double * const cw2 = this->up_w2_.get();

            this->up_kerSummator_->nullify_Ieta_2();
            for (size_t n2_psi = 0; n2_psi < this->N2(); ++n2_psi) {

                const double Psi_2 = psi_b_mnus_a * cz2[n2_psi] + psi_a_plus_b;
                update_psi_trigonometry_(Psi_2);
                calc_I_lam_N3_nfuncs_wrapper_(cos_psi_2_ * cw2[n2_psi]);
            } // N2 for
            this->up_kerSummator_->multiply_Ieta_2(psi_b_mnus_a);
            this->up_kerSummator_->accumulate_Ipsi_1(cw1[n_theta_1]);
        } // N1 for
        this->up_kerSummator_->multiply_Ipsi_1(theta_b_mnus_a);
        this->up_kerSummator_->assign_Ipsi_1_to( &Isub_[this->up_kernel_->size() * m]);
    }  // 7 or 9 == m for triangular or quadrilateral elements
    gather_Iss_();
}

template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::gather_Iss_() noexcept {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Collect sub ranges for each kernel
    this->Iss_tot_ = dcomplex(0.0, 0.0);
    for (size_t i = 0; i < t_ker_sz; ++i) {
        dcomplex Iss_isum(0.0, 0.0);
        for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {
            Iss_isum += Isub_[i + t_ker_sz * m];
        }
        this->Iss_[i] = Iss_isum * this->up_kernel_->precomputed_jacobian();
        this->Iss_tot_ += this->Iss_[i];
    }
}

template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::update_theta_trigonometry_(const double Theta_1) noexcept {

    cos_theta_1_ = cos(Theta_1);
    sin_theta_1_ = sin(Theta_1);
    tan_theta_1_ = tan(Theta_1);
    tan_05pi_mns_theta_1_ = tan(0.5 * M_PI - Theta_1);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::update_psi_trigonometry_(const double Psi_2) noexcept {

    cos_psi_2_ = cos(Psi_2);
    sin_psi_2_ = sin(Psi_2);
    tan_psi_2_ = tan(Psi_2);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::single_n_func_I_(const double mlt_factor) noexcept {

    setup_limits34_pntr_to_single_();
    calc_n_func_();
    // No changes in the algorithm
    this->up_kerSummator_->accumulate_Ieta_2(mlt_factor);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::double_n_func_II_(const double mlt_factor) noexcept {

    setup_limits34_pntr_to_double_first_();
    calc_n_func_();
    // Have to consider two terms: the first
    this->up_kerSummator_->accumulate_Ieta_2(mlt_factor);

    setup_limits34_pntr_to_double_secnd_();
    calc_n_func_();
    // Have to consider two terms: the second
    this->up_kerSummator_->accumulate_Ieta_2(mlt_factor);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::calc_n_func_() noexcept {

    double ueta_a = 0.0, ueta_b = 0.0;
    get_ueta_limits3_nfun_(ueta_a, ueta_b);
    const double eta_plus_ba = 0.5 * (ueta_b + ueta_a);
    const double eta_mnus_ba = 0.5 * (ueta_b - ueta_a);

    const double * const cw3 = this->up_w3_.get();
    const double * const cz3 = this->up_z3_.get();

    this->up_kerSummator_->nullify_Ilam_3();
    for (size_t n3_eta = 0; n3_eta < this->N3(); ++n3_eta) {

        const double Eta_3 = eta_mnus_ba * cz3[n3_eta] + eta_plus_ba;
        calc_afuncs_(Eta_3);
        this->up_kerSummator_->accumulate_Ilam_3(cw3[n3_eta]);
    }
    this->up_kerSummator_->multiply_Ilam_3(eta_mnus_ba);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_EA<ParticularKernel>::calc_afuncs_(const double Eta_3) noexcept {

    const double lam_lim = get_lam_limit_(Eta_3);
    const double J_lam = 0.5 * lam_lim;

    const double * const cz4 = this->up_z4_.get();
    const double * const cw4 = this->up_w4_.get();

    this->up_kerSummator_->nullify_Irho_4();
    for (size_t n4_lam = 0; n4_lam < this->N4(); ++n4_lam) {

        const double Lam_4 = J_lam * (cz4[n4_lam] + 1.0);
        double uvxi_p[3], uvxi_q[3];
        ea_calc_pq_simplex_(uvxi_p, uvxi_q, Eta_3, Lam_4); // out, out, in

        this->up_kernel_->update_rp(uvxi_p);
        this->up_kernel_->update_rq(uvxi_q);
        // Precaches the Greens function here.
        this->up_kernel_->precompute_rp_rq_data();
        this->up_kerSummator_->accumulate_Irho_4(cw4[n4_lam] * Lam_4 * Lam_4);
    }
    this->up_kerSummator_->multiply_Irho_4(J_lam);
}

// Instantiation of the Triangular Constant Kernels
template class DirectfnAlgorithm_EA<TriangularKernel_Constant_ST>;
template class DirectfnAlgorithm_EA<TriangularKernel_Constant_EA>;
template class DirectfnAlgorithm_EA<TriangularKernel_Constant_VA>;

template class DirectfnAlgorithm_EA<TriangularKernel_RWG_WS>;
template class DirectfnAlgorithm_EA<TriangularKernel_RWG_SS>;
template class DirectfnAlgorithm_EA<TriangularKernel_nxRWG_SS>;

// Instantiation of the Quadrilateral Planar Kernels
template class DirectfnAlgorithm_EA<QuadrilateralKernel_PlanarScalar>;
template class DirectfnAlgorithm_EA<QuadrilateralKernel_PlanarVectorWS>;
template class DirectfnAlgorithm_EA<QuadrilateralKernel_PlanarVectorSS>;

template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class DirectfnAlgorithm_EA<QuadKer_PlanVH_VolKer3_KerTyp4>;

// Instantiation of the Quadrilateral Curvilinear Kernels
template class DirectfnAlgorithm_EA<QuadrilateralKernel_CurvilinearScalar>;
template class DirectfnAlgorithm_EA<QuadrilateralKernel_CurvilinearVectorWS>;
template class DirectfnAlgorithm_EA<QuadrilateralKernel_CurvilinearVectorSS>;

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Triangular_EA<ParticularKernel>::Triangular_EA():
DirectfnAlgorithm_EA<ParticularKernel>(),
pf_thetalim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_psilim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
fptr_nfunc_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
eta_pq_sgn_(0.0),
pf_eta3_lim_crnt_(nullptr),
pf_eta3_lim_single_crnt_(nullptr),
pf_eta3_lim_double_frst_crnt_(nullptr),
pf_eta3_lim_double_scnd_crnt_(nullptr),
pf_eta3_lim_single_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_eta3_lim_double_frst_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_eta3_lim_double_scnd_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_lam4_lim_crnt_(nullptr),
pf_lam4_lim_single_crnt_(nullptr),
pf_lam4_lim_double_frst_crnt_(nullptr),
pf_lam4_lim_double_scnd_crnt_(nullptr),
pf_lam4_lim_single_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_lam4_lim_double_frst_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_lam4_lim_double_scnd_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr} {

    allocate_I_vars_();
    initialize_limit_fptrs_();
}


//virtual
template <typename ParticularKernel>
Triangular_EA<ParticularKernel>::~Triangular_EA() {
}

template <typename ParticularKernel>
string Triangular_EA<ParticularKernel>::name() const noexcept {
    return string("Triangular_EA_Constant algorithm");
}

template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::allocate_I_vars_() noexcept {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_ranges_6_() * t_ker_sz]);
}

template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::initialize_limit_fptrs_() noexcept {
    // N1
    pf_thetalim_arr_[0] = triag_ea_theta_limits<0>;
    pf_thetalim_arr_[1] = triag_ea_theta_limits<1>;
    pf_thetalim_arr_[2] = triag_ea_theta_limits<2>;
    pf_thetalim_arr_[3] = triag_ea_theta_limits<3>;
    pf_thetalim_arr_[4] = triag_ea_theta_limits<4>;
    pf_thetalim_arr_[5] = triag_ea_theta_limits<5>;
    // N2
    pf_psilim_arr_[0] = triag_ea_psi_limits<0>;
    pf_psilim_arr_[1] = triag_ea_psi_limits<1>;
    pf_psilim_arr_[2] = triag_ea_psi_limits<2>;
    pf_psilim_arr_[3] = triag_ea_psi_limits<3>;
    pf_psilim_arr_[4] = triag_ea_psi_limits<4>;
    pf_psilim_arr_[5] = triag_ea_psi_limits<5>;
    // N3
    fptr_nfunc_arr_[0] = &DirectfnAlgorithm_EA<ParticularKernel>::double_n_func_II_;
    fptr_nfunc_arr_[1] = &DirectfnAlgorithm_EA<ParticularKernel>::double_n_func_II_;
    fptr_nfunc_arr_[2] = &DirectfnAlgorithm_EA<ParticularKernel>::double_n_func_II_;
    fptr_nfunc_arr_[3] = &DirectfnAlgorithm_EA<ParticularKernel>::single_n_func_I_;
    fptr_nfunc_arr_[4] = &DirectfnAlgorithm_EA<ParticularKernel>::single_n_func_I_;
    fptr_nfunc_arr_[5] = &DirectfnAlgorithm_EA<ParticularKernel>::single_n_func_I_;

    // Eta 3 : single
    pf_eta3_lim_single_arr_[0] = nullptr;
    pf_eta3_lim_single_arr_[1] = nullptr;
    pf_eta3_lim_single_arr_[2] = nullptr;
    pf_eta3_lim_single_arr_[3] = &Triangular_EA::eta3_lims_single_;
    pf_eta3_lim_single_arr_[4] = &Triangular_EA::eta3_lims_single_;
    pf_eta3_lim_single_arr_[5] = &Triangular_EA::eta3_lims_single_;
    // Eta 3 : double - I (upper)
    pf_eta3_lim_double_frst_arr_[0] = &Triangular_EA::eta3_lims_double_frst_I_12;
    pf_eta3_lim_double_frst_arr_[1] = &Triangular_EA::eta3_lims_double_frst_I_12;
    pf_eta3_lim_double_frst_arr_[2] = &Triangular_EA::eta3_lims_double_frst_I_3;
    pf_eta3_lim_double_frst_arr_[3] = nullptr;
    pf_eta3_lim_double_frst_arr_[4] = nullptr;
    pf_eta3_lim_double_frst_arr_[5] = nullptr;
    // Eta 3 : double - II (second)
    pf_eta3_lim_double_scnd_arr_[0] = &Triangular_EA::eta3_lims_double_scnd_II_12_;
    pf_eta3_lim_double_scnd_arr_[1] = &Triangular_EA::eta3_lims_double_scnd_II_12_;
    pf_eta3_lim_double_scnd_arr_[2] = &Triangular_EA::eta3_lims_double_scnd_II_3_;
    pf_eta3_lim_double_scnd_arr_[3] = nullptr;
    pf_eta3_lim_double_scnd_arr_[4] = nullptr;
    pf_eta3_lim_double_scnd_arr_[5] = nullptr;

    // Lam 4 : single
    pf_lam4_lim_single_arr_[0] = nullptr;
    pf_lam4_lim_single_arr_[1] = nullptr;
    pf_lam4_lim_single_arr_[2] = nullptr;
    pf_lam4_lim_single_arr_[3] = &Triangular_EA::lambda_L1_;
    pf_lam4_lim_single_arr_[4] = &Triangular_EA::lambda_L2_;
    pf_lam4_lim_single_arr_[5] = &Triangular_EA::lambda_L2_;
    // Lam 4 : double - I
    pf_lam4_lim_double_frst_arr_[0] = &Triangular_EA::lambda_L3_;
    pf_lam4_lim_double_frst_arr_[1] = &Triangular_EA::lambda_L3_;
    pf_lam4_lim_double_frst_arr_[2] = &Triangular_EA::lambda_L3_;
    pf_lam4_lim_double_frst_arr_[3] = nullptr;
    pf_lam4_lim_double_frst_arr_[4] = nullptr;
    pf_lam4_lim_double_frst_arr_[5] = nullptr;
    // Lam 4 : double - II
    pf_lam4_lim_double_scnd_arr_[0] = &Triangular_EA::lambda_L2_abs_;
    pf_lam4_lim_double_scnd_arr_[1] = &Triangular_EA::lambda_L2_abs_;
    pf_lam4_lim_double_scnd_arr_[2] = &Triangular_EA::lambda_L1_;
    pf_lam4_lim_double_scnd_arr_[3] = nullptr;
    pf_lam4_lim_double_scnd_arr_[4] = nullptr;
    pf_lam4_lim_double_scnd_arr_[5] = nullptr;
}

template <typename ParticularKernel>
void  Triangular_EA<ParticularKernel>::eta3_lims_single_(double & Ua, double & Ub) noexcept {
    triag_ea_eta_limits(Ua, Ub);
}

template <typename ParticularKernel>
void  Triangular_EA<ParticularKernel>::eta3_lims_double_frst_I_12(double & Ua, double & Ub) noexcept {

    const double  gamma = gamma_sin_cos_();
    const double  Psi_eta = (1.0 - gamma) / (1.0 + gamma);

    Ua = 0;
    Ub = Psi_eta;
}

template <typename ParticularKernel>
void  Triangular_EA<ParticularKernel>::eta3_lims_double_frst_I_3(double & Ua, double & Ub) noexcept {

    Ua = 0.0;
    Ub = gamma_tan_();
}

template <typename ParticularKernel>
void  Triangular_EA<ParticularKernel>::eta3_lims_double_scnd_II_12_(double & Ua, double & Ub) noexcept {

    const double  gamma = gamma_sin_cos_();
    const double  Psi_eta = (1.0 - gamma) / (1.0 + gamma);

    Ua = Psi_eta;
    Ub = 1.0;
}

template <typename ParticularKernel>
void  Triangular_EA<ParticularKernel>::eta3_lims_double_scnd_II_3_(double & Ua, double & Ub) noexcept {

    Ua = gamma_tan_();
    Ub = 1.0;
}

template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::gamma_sin_cos_() const noexcept {

    return this->tan_psi_2_ / (this->sin_theta_1_ + sqrt(3.0) * this->cos_theta_1_);
}

template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::gamma_tan_() const noexcept {

    return  sqrt(3.0) / this->tan_theta_1_;
}

template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::lambda_L1_(const double Eta_3) const noexcept {

    const double sqrt3 = sqrt(3.0);
    const double ddenom = this->cos_psi_2_ * (this->sin_theta_1_ / sqrt3 - this->cos_theta_1_);
    return (1.0 - fabs(Eta_3)) /ddenom;
}

template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::lambda_L2_(const double Eta_3) const noexcept {

    return  sqrt(3.0) * (1.0 - Eta_3) / this->sin_psi_2_;
}

template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::lambda_L2_abs_(const double Eta_3) const noexcept {

    return  sqrt(3.0) * (1.0 - fabs(Eta_3)) / this->sin_psi_2_;
}

template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::lambda_L3_(const double Eta_3) const noexcept {

    const double sqrt3 = sqrt(3.0);
    const double ddenom = this->cos_psi_2_ * (this->sin_theta_1_ / sqrt3 + this->cos_theta_1_);
    return (1.0 + Eta_3) /ddenom;
}

//virtual
template <typename ParticularKernel>
size_t Triangular_EA<ParticularKernel>::sub_ranges_numb_m_() const noexcept {
    return sub_ranges_6_();
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::update_theta_lims_N1_fptr_(const size_t m) noexcept {
    this->pf_thetalim_crnt_ = pf_thetalim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::update_psi_lims_N2_fptr_(const size_t m) noexcept {
    this->pf_psilim_crnt_ = pf_psilim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::update_eta_lims_N3_fptr_(const size_t m) noexcept {

    this->fptr_nfunc_crnt_ = fptr_nfunc_arr_[m];
    // Single    case
    this->pf_eta3_lim_single_crnt_ = pf_eta3_lim_single_arr_[m];
    // Double I  case
    this->pf_eta3_lim_double_frst_crnt_ = pf_eta3_lim_double_frst_arr_[m];
    // Double II case
    this->pf_eta3_lim_double_scnd_crnt_ = pf_eta3_lim_double_scnd_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::update_lam_lims_N4_fprt_(const size_t m) noexcept {

    // Single    case
    this->pf_lam4_lim_single_crnt_      = pf_lam4_lim_single_arr_[m];
    // Double I  case
    this->pf_lam4_lim_double_frst_crnt_ = pf_lam4_lim_double_frst_arr_[m];
    // Double II case
    this->pf_lam4_lim_double_scnd_crnt_ = pf_lam4_lim_double_scnd_arr_[m];
}

//virtual
template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::init_psi_a_by_theta_() noexcept {
    return atan(this->sin_theta_1_ - sqrt(3.0) * this->cos_theta_1_);
}

//virtual
template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::init_psi_b_by_theta_() noexcept {
    return atan(this->sin_theta_1_ + sqrt(3.0) * this->cos_theta_1_);
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::calc_I_lam_N3_nfuncs_wrapper_(const double mlt_factor) noexcept {

    // This parameter is used for triangular simplex only
    eta_pq_sgn_ = 1.0;
    // n_function pointer can be pointer to a double call of a_function
    // or to a single call of a_function
    (this->*this->fptr_nfunc_crnt_)(mlt_factor);

    // This parameter is used for triangular simplex only
    eta_pq_sgn_ = -1.0;
    (this->*this->fptr_nfunc_crnt_)(mlt_factor);
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::setup_limits34_pntr_to_single_() noexcept {

    this->pf_eta3_lim_crnt_ = pf_eta3_lim_single_crnt_;
    this->pf_lam4_lim_crnt_ = pf_lam4_lim_single_crnt_;
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::setup_limits34_pntr_to_double_first_() noexcept {

    pf_eta3_lim_crnt_ = pf_eta3_lim_double_frst_crnt_;
    pf_lam4_lim_crnt_ = pf_lam4_lim_double_frst_crnt_;
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::setup_limits34_pntr_to_double_secnd_() noexcept {

    pf_eta3_lim_crnt_ = pf_eta3_lim_double_scnd_crnt_;
    pf_lam4_lim_crnt_ = pf_lam4_lim_double_scnd_crnt_;
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::get_ueta_limits3_nfun_(double & Ua, double & Ub) noexcept {

    (this->*pf_eta3_lim_crnt_)(Ua, Ub);
}

//virtual
template <typename ParticularKernel>
double Triangular_EA<ParticularKernel>::get_lam_limit_(const double Eta_3) noexcept {

    return (this->*pf_lam4_lim_crnt_)(Eta_3);
}

//virtual
template <typename ParticularKernel>
void Triangular_EA<ParticularKernel>::ea_calc_pq_simplex_(double xi_p_out[3], double xi_q_out[3],
                                                const double Eta3, const double Lam4) noexcept {
    const double eta_p = Eta3 * this->eta_pq_sgn_;
    const double xi_p  = Lam4 * this->sin_psi_2_;
    triag_st_make_simplex(xi_p_out, eta_p, xi_p);

    const double eta_q = this->eta_pq_sgn_ * (Lam4 * this->cos_psi_2_ * this->cos_theta_1_ - Eta3);
    const double xi_q  = Lam4 * this->cos_psi_2_ * this->sin_theta_1_;
    triag_st_make_simplex(xi_q_out, eta_q, xi_q);
}


template class Triangular_EA<TriangularKernel_Constant_ST>;
template class Triangular_EA<TriangularKernel_Constant_EA>;
template class Triangular_EA<TriangularKernel_Constant_VA>;

template class Triangular_EA<TriangularKernel_RWG_WS>;
template class Triangular_EA<TriangularKernel_RWG_SS>;
template class Triangular_EA<TriangularKernel_nxRWG_SS>;

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Quadrilateral_EA<ParticularKernel>::Quadrilateral_EA() :
DirectfnAlgorithm_EA<ParticularKernel>(),
pf_thetalim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_psilim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
fptr_nfunc_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_eta3_lim_crnt_(nullptr),
pf_eta3_lim_single_crnt_(nullptr),
pf_eta3_lim_double_frst_crnt_(nullptr),
pf_eta3_lim_double_scnd_crnt_(nullptr),
pf_eta3_lim_single_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_eta3_lim_double_frst_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_eta3_lim_double_scnd_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_lam4_lim_crnt_(nullptr),
pf_lam4_lim_single_crnt_(nullptr),
pf_lam4_lim_double_frst_crnt_(nullptr),
pf_lam4_lim_double_scnd_crnt_(nullptr),
pf_lam4_lim_single_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_lam4_lim_double_frst_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_lam4_lim_double_scnd_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr} {

    // The kernel pointer as well as the kernel-array-summator pointer
    // has been allocated in the DirectfnInterface<ParticularKernel>

    // Here the memory for arrays is allocated according to kernel size.
    allocate_I_vars_();
    initialize_limit_fptrs_();
}

//virtual
template <typename ParticularKernel>
Quadrilateral_EA<ParticularKernel>::~Quadrilateral_EA() {
}

//virtual
template <typename ParticularKernel>
string Quadrilateral_EA<ParticularKernel>::name() const noexcept {
    return string("Quadrilateral_EA");
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::allocate_I_vars_() noexcept {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_ranges_8_() * t_ker_sz]);
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::initialize_limit_fptrs_() noexcept {
    // N1
    pf_thetalim_arr_[0] = quad_ea_theta_limits<0>;
    pf_thetalim_arr_[1] = quad_ea_theta_limits<1>;
    pf_thetalim_arr_[2] = quad_ea_theta_limits<2>;
    pf_thetalim_arr_[3] = quad_ea_theta_limits<3>;
    pf_thetalim_arr_[4] = quad_ea_theta_limits<4>;
    pf_thetalim_arr_[5] = quad_ea_theta_limits<5>;
    pf_thetalim_arr_[6] = quad_ea_theta_limits<6>;
    pf_thetalim_arr_[7] = quad_ea_theta_limits<7>;
    // N2
    pf_psilim_arr_[0] = quad_ea_psi_limits<0>;
    pf_psilim_arr_[1] = quad_ea_psi_limits<1>;
    pf_psilim_arr_[2] = quad_ea_psi_limits<2>;
    pf_psilim_arr_[3] = quad_ea_psi_limits<3>;
    pf_psilim_arr_[4] = quad_ea_psi_limits<4>;
    pf_psilim_arr_[5] = quad_ea_psi_limits<5>;
    pf_psilim_arr_[6] = quad_ea_psi_limits<6>;
    pf_psilim_arr_[7] = quad_ea_psi_limits<7>;
    // N3
    //fptr_nfunc_arr_[0] = &DirectfnAlgorithm_EA::single_n_func_I_;
    fptr_nfunc_arr_[0] = &Quadrilateral_EA::single_n_func_I_;
    fptr_nfunc_arr_[1] = &Quadrilateral_EA::double_n_func_II_;
    fptr_nfunc_arr_[2] = &Quadrilateral_EA::double_n_func_II_;
    fptr_nfunc_arr_[3] = &Quadrilateral_EA::double_n_func_II_;
    fptr_nfunc_arr_[4] = &Quadrilateral_EA::double_n_func_II_;
    fptr_nfunc_arr_[5] = &Quadrilateral_EA::double_n_func_II_;
    fptr_nfunc_arr_[6] = &Quadrilateral_EA::single_n_func_I_;
    fptr_nfunc_arr_[7] = &Quadrilateral_EA::double_n_func_II_;

    // Eta 3 : single
    pf_eta3_lim_single_arr_[0] = &Quadrilateral_EA::eta3_lims_single_;
    pf_eta3_lim_single_arr_[1] = nullptr;
    pf_eta3_lim_single_arr_[2] = nullptr;
    pf_eta3_lim_single_arr_[3] = nullptr;
    pf_eta3_lim_single_arr_[4] = nullptr;
    pf_eta3_lim_single_arr_[5] = nullptr;
    pf_eta3_lim_single_arr_[6] = &Quadrilateral_EA::eta3_lims_single_;
    pf_eta3_lim_single_arr_[7] = nullptr;
    // Eta 3 : double - I (upper)
    pf_eta3_lim_double_frst_arr_[0] = nullptr;
    pf_eta3_lim_double_frst_arr_[1] = &Quadrilateral_EA::eta3_lims_double_m1Fm_;
    pf_eta3_lim_double_frst_arr_[2] = &Quadrilateral_EA::eta3_lims_double_m1Tm_;
    pf_eta3_lim_double_frst_arr_[3] = &Quadrilateral_EA::eta3_lims_double_m1Fm_;
    pf_eta3_lim_double_frst_arr_[4] = &Quadrilateral_EA::eta3_lims_double_m1Tp_;
    pf_eta3_lim_double_frst_arr_[5] = &Quadrilateral_EA::eta3_lims_double_m1Fp_;
    pf_eta3_lim_double_frst_arr_[6] = nullptr;
    pf_eta3_lim_double_frst_arr_[7] = &Quadrilateral_EA::eta3_lims_double_Fpp1_;
    // Eta 3 : double - II (second)
    pf_eta3_lim_double_scnd_arr_[0] = nullptr;
    pf_eta3_lim_double_scnd_arr_[1] = &Quadrilateral_EA::eta3_lims_double_Fmp1_;
    pf_eta3_lim_double_scnd_arr_[2] = &Quadrilateral_EA::eta3_lims_double_Tmp1_;
    pf_eta3_lim_double_scnd_arr_[3] = &Quadrilateral_EA::eta3_lims_double_Fmp1_;
    pf_eta3_lim_double_scnd_arr_[4] = &Quadrilateral_EA::eta3_lims_double_Tpp1_;
    pf_eta3_lim_double_scnd_arr_[5] = &Quadrilateral_EA::eta3_lims_double_Fpp1_;
    pf_eta3_lim_double_scnd_arr_[6] = nullptr;
    pf_eta3_lim_double_scnd_arr_[7] = &Quadrilateral_EA::eta3_lims_double_m1Fp_;

    // Lam 4 : single
    pf_lam4_lim_single_arr_[0] = &Quadrilateral_EA::Lam_xp1_cc_;
    pf_lam4_lim_single_arr_[1] = nullptr;
    pf_lam4_lim_single_arr_[2] = nullptr;
    pf_lam4_lim_single_arr_[3] = nullptr;
    pf_lam4_lim_single_arr_[4] = nullptr;
    pf_lam4_lim_single_arr_[5] = nullptr;
    pf_lam4_lim_single_arr_[6] = &Quadrilateral_EA::Lam_xm1_cc_;
    pf_lam4_lim_single_arr_[7] = nullptr;
    // Lam 4 : double - I upper first
    pf_lam4_lim_double_frst_arr_[0] = nullptr;
    pf_lam4_lim_double_frst_arr_[1] = &Quadrilateral_EA::Lam_xp1_cc_;
    pf_lam4_lim_double_frst_arr_[2] = &Quadrilateral_EA::Lam_xp1_cc_;
    pf_lam4_lim_double_frst_arr_[3] = &Quadrilateral_EA::Lam_xp1_cc_;
    pf_lam4_lim_double_frst_arr_[4] = &Quadrilateral_EA::Lam_0p2_sc_;
    pf_lam4_lim_double_frst_arr_[5] = &Quadrilateral_EA::Lam_0p2_0s_;
    pf_lam4_lim_double_frst_arr_[6] = nullptr;
    pf_lam4_lim_double_frst_arr_[7] = &Quadrilateral_EA::Lam_xm1_cc_;
    // Lam 4 : double - II second
    pf_lam4_lim_double_scnd_arr_[0] = nullptr;
    pf_lam4_lim_double_scnd_arr_[1] = &Quadrilateral_EA::Lam_0p2_0s_;
    pf_lam4_lim_double_scnd_arr_[2] = &Quadrilateral_EA::Lam_0p2_sc_;
    pf_lam4_lim_double_scnd_arr_[3] = &Quadrilateral_EA::Lam_0p2_0s_;
    pf_lam4_lim_double_scnd_arr_[4] = &Quadrilateral_EA::Lam_xm1_cc_;
    pf_lam4_lim_double_scnd_arr_[5] = &Quadrilateral_EA::Lam_xm1_cc_;
    pf_lam4_lim_double_scnd_arr_[6] = nullptr;
    pf_lam4_lim_double_scnd_arr_[7] = &Quadrilateral_EA::Lam_0p2_0s_;
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_single_(double & Ua, double & Ub) noexcept {
    quadr_ea_eta_limits(Ua, Ub);
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::F_plus_1_() const noexcept {
    return 2.0 * this->cos_theta_1_ / this->tan_psi_2_ + 1.0;
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::F_mnus_1_() const noexcept {
    return 2.0 * this->cos_theta_1_ / this->tan_psi_2_ - 1.0;
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::T_plus_1_() const noexcept {
    return 2.0 * this->tan_05pi_mns_theta_1_ + 1.0;
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::T_mnus_1_() const noexcept {
    return 2.0 * this->tan_05pi_mns_theta_1_ - 1.0;
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::Lam_xp1_cc_(const double x) const noexcept {
    return (x + 1.0) / (this->cos_theta_1_ * this->cos_psi_2_);
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::Lam_xm1_cc_(const double x) const noexcept {
    return (x - 1.0) / (this->cos_theta_1_ * this->cos_psi_2_);
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::Lam_0p2_sc_(const double ) const noexcept {
    return 2.0 / (this->sin_theta_1_ * this->cos_psi_2_);
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::Lam_0p2_0s_(const double ) const noexcept {
    return 2.0 / this->sin_psi_2_;
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_m1Fm_(double & Ua, double & Ub) noexcept {
    Ua = -1.0;
    Ub =  F_mnus_1_();
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_m1Tm_(double & Ua, double & Ub) noexcept {
    Ua = -1.0;
    Ub =  T_mnus_1_();
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_m1Fp_(double & Ua, double & Ub) noexcept {
    Ua = -1.0;
    Ub =  F_plus_1_();
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_m1Tp_(double & Ua, double & Ub) noexcept {
    Ua = -1.0;
    Ub =  T_plus_1_();
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_Fmp1_(double & Ua, double & Ub) noexcept {
    Ua = F_mnus_1_();
    Ub = 1.0;
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_Tmp1_(double & Ua, double & Ub) noexcept {
    Ua = T_mnus_1_();
    Ub = 1.0;
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_Fpp1_(double & Ua, double & Ub) noexcept {
    Ua = F_plus_1_();
    Ub = 1.0;
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::eta3_lims_double_Tpp1_(double & Ua, double & Ub) noexcept {
    Ua = T_plus_1_();
    Ub = 1.0;
}

//virtual
template <typename ParticularKernel>
size_t Quadrilateral_EA<ParticularKernel>::sub_ranges_numb_m_() const noexcept {
    return  sub_ranges_8_();
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::update_theta_lims_N1_fptr_(const size_t m) noexcept {
    this->pf_thetalim_crnt_ = this->pf_thetalim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::update_psi_lims_N2_fptr_(const size_t m) noexcept {
    this->pf_psilim_crnt_ = this->pf_psilim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::update_eta_lims_N3_fptr_(const size_t m) noexcept {

    this->fptr_nfunc_crnt_ = this->fptr_nfunc_arr_[m];
    // Single    case
    this->pf_eta3_lim_single_crnt_ = this->pf_eta3_lim_single_arr_[m];
    // Double I  case
    this->pf_eta3_lim_double_frst_crnt_ = this->pf_eta3_lim_double_frst_arr_[m];
    // Double II case
    this->pf_eta3_lim_double_scnd_crnt_ = this->pf_eta3_lim_double_scnd_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::update_lam_lims_N4_fprt_(const size_t m) noexcept {

    // Single    case
    this->pf_lam4_lim_single_crnt_      = this->pf_lam4_lim_single_arr_[m];
    // Double I  case
    this->pf_lam4_lim_double_frst_crnt_ = this->pf_lam4_lim_double_frst_arr_[m];
    // Double II case
    this->pf_lam4_lim_double_scnd_crnt_ = this->pf_lam4_lim_double_scnd_arr_[m];
}

//virtual
template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::init_psi_a_by_theta_() noexcept {
    return atan(this->cos_theta_1_);
}

//virtual
template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::init_psi_b_by_theta_() noexcept {
    return atan(this->sin_theta_1_);
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::calc_I_lam_N3_nfuncs_wrapper_(const double mlt_factor) noexcept {
    // n_function pointer can be pointer to a double call of a_function
    // or to a single call of a_function
    (this->*this->fptr_nfunc_crnt_)(mlt_factor);
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::setup_limits34_pntr_to_single_() noexcept {

    pf_eta3_lim_crnt_ = pf_eta3_lim_single_crnt_;
    pf_lam4_lim_crnt_ = pf_lam4_lim_single_crnt_;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::setup_limits34_pntr_to_double_first_() noexcept {

    pf_eta3_lim_crnt_ = pf_eta3_lim_double_frst_crnt_;
    pf_lam4_lim_crnt_ = pf_lam4_lim_double_frst_crnt_;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::setup_limits34_pntr_to_double_secnd_() noexcept {

    pf_eta3_lim_crnt_ = pf_eta3_lim_double_scnd_crnt_;
    pf_lam4_lim_crnt_ = pf_lam4_lim_double_scnd_crnt_;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::get_ueta_limits3_nfun_(double & Ua, double & Ub) noexcept {
    (this->*pf_eta3_lim_crnt_)(Ua, Ub);
}

template <typename ParticularKernel>
double Quadrilateral_EA<ParticularKernel>::get_lam_limit_(const double Eta_3) noexcept {
    return (this->*pf_lam4_lim_crnt_)(Eta_3);
}

template <typename ParticularKernel>
void Quadrilateral_EA<ParticularKernel>::ea_calc_pq_simplex_(double uvxi_p[3], double uvxi_q[3],
                                                const double Eta_3, const double Lam_4) noexcept {
    uvxi_p[0] = Eta_3;
    uvxi_p[1] = Lam_4 * this->sin_psi_2_ - 1.0;

    uvxi_q[0] = Lam_4 * this->cos_theta_1_ * this->cos_psi_2_ - Eta_3;
    uvxi_q[1] = Lam_4 * this->sin_theta_1_ * this->cos_psi_2_ - 1.0;
}


// Instantiation of the Quadrilateral Constant Kernels
template class Quadrilateral_EA<QuadrilateralKernel_PlanarScalar>;
template class Quadrilateral_EA<QuadrilateralKernel_PlanarVectorWS>;
template class Quadrilateral_EA<QuadrilateralKernel_PlanarVectorSS>;

template class Quadrilateral_EA<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class Quadrilateral_EA<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class Quadrilateral_EA<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class Quadrilateral_EA<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class Quadrilateral_EA<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class Quadrilateral_EA<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class Quadrilateral_EA<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class Quadrilateral_EA<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class Quadrilateral_EA<QuadrilateralKernel_CurvilinearScalar>;
template class Quadrilateral_EA<QuadrilateralKernel_CurvilinearVectorWS>;
template class Quadrilateral_EA<QuadrilateralKernel_CurvilinearVectorSS>;

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file

