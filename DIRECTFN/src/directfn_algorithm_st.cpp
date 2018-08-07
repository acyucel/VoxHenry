#include <iostream>
#include <iomanip>
#include "directfn_algorithm_st.h"
#include "directfn_kernel_tri.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"
#include "directfn_kernel_quad_voxhenry.h"

using  std::cout;
using  std::endl;
using  std::setprecision;
using  std::setw;
using  std::to_string;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
DirectfnAlgorithm_ST<ParticularKernel>::DirectfnAlgorithm_ST():
DirectfnInterface<ParticularKernel>(),
sin_PSI_(0.0),
cos_PSI_(0.0),
t_Eta_2_(0.0),
Isub_(nullptr),
Isd_k_(nullptr),
pf_subtrg_crnt_(nullptr),
pf_psilim_crnt_(nullptr),
pf_etalim_crnt_(nullptr) {

}

//virtual
template <typename ParticularKernel>
DirectfnAlgorithm_ST<ParticularKernel>::~DirectfnAlgorithm_ST() {
}

//virtual
template <typename ParticularKernel>
void  DirectfnAlgorithm_ST<ParticularKernel>::precompute_psi_trigonomety_(const double t_Psi) noexcept {

    sin_PSI_ = sin(t_Psi);
    cos_PSI_ = cos(t_Psi);
}

//virtual
template <typename ParticularKernel>
void DirectfnAlgorithm_ST<ParticularKernel>::do_I_surface_surface_() {

    // Preliminarily routine
    init_Isub_();
    // Go throw three triangular / four quadrilateral sub-elements deviding the whole domain
    for (size_t kk = 0; kk < sub_figures_numb_k_(); ++kk) {

        update_subtri_fptrs_(kk);
        // Go throw subranges
        for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {

            // Update limits f-pointers for a given subrange
            update_psilims_fptrs_N1_(m);
            update_etalims_fptrs_N2_(m);
            update_lamlims_fptrs_N3_(m);

            double Psi_A = 0.0, Psi_B = 0.0;
            pf_psilim_crnt_(Psi_A, Psi_B);
            const double half_psi_a_plus_b = 0.5 * (Psi_A + Psi_B);
            const double half_psi_b_mnus_a = 0.5 * (Psi_B - Psi_A);

            const double * const cz1 = this->up_z1_.get();
            const double * const cw1 = this->up_w1_.get();

            this->up_kerSummator_->nullify_Ipsi_1();
            for (size_t n1_psi = 0; n1_psi < this->N1(); ++n1_psi) {

                const double t_Psi = half_psi_b_mnus_a * cz1[n1_psi] + half_psi_a_plus_b;
                precompute_psi_trigonomety_(t_Psi);

                double Eta_a = 0.0, Eta_b = 0.0;
                pf_etalim_crnt_(t_Psi, Eta_a, Eta_b);
                const double half_eta_b_mnus_a = 0.5 * (Eta_b - Eta_a);
                const double half_eta_a_plus_b = 0.5 * (Eta_a + Eta_b);

                const double * const cz2 = this->up_z2_.get();
                const double * const cw2 = this->up_w2_.get();

                this->up_kerSummator_->nullify_Ieta_2();
                for (size_t n2_eta = 0; n2_eta < this->N2(); ++n2_eta) {

                    t_Eta_2_ = half_eta_b_mnus_a * cz2[n2_eta] + half_eta_a_plus_b;
                    calc_n_functions_();
                    this->up_kerSummator_->accumulate_Ieta_2(cw2[n2_eta]);
                }  // N2 for
                this->up_kerSummator_->multiply_Ieta_2(half_eta_b_mnus_a * sin_PSI_);
                this->up_kerSummator_->accumulate_Ipsi_1(cw1[n1_psi]);
            } // N1 for
            this->up_kerSummator_->multiply_Ipsi_1(half_psi_b_mnus_a);
            this->up_kerSummator_->assign_Ipsi_1_to( &Isd_k_[this->up_kernel_->size() * m] );
        } // 8 or 6  ==  m for                          // Isd_k_[i + ker_sz * m] = I_psi[i];
        gather_Isub_(kk);
    } // 3 or 4 == kk for triangular or quadrilateral elements
    gather_Iss_();  // Iss_[i] and Iss_tot
}

template <typename ParticularKernel>
void DirectfnAlgorithm_ST<ParticularKernel>::init_Isub_() noexcept {

    const size_t t_ker_sz = this->up_kernel_->size();

    // Clear tmp buffer for subdomain integrals for quadrilateral
    for (size_t kk = 0; kk < sub_figures_numb_k_(); ++kk) {
        for (size_t i = 0; i < t_ker_sz; ++i) {
            Isub_[i + t_ker_sz * kk] = dcomplex(0., 0.);
        }
    }
}

template <typename ParticularKernel>
void DirectfnAlgorithm_ST<ParticularKernel>::gather_Isub_(const size_t kk) noexcept {

    const size_t ker_sz = this->up_kernel_->size();

    for (size_t i = 0; i < ker_sz; ++i) {
        for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {
            Isub_[i + ker_sz * kk] += Isd_k_[i + ker_sz * m];
        }
    }
}

template <typename ParticularKernel>
void DirectfnAlgorithm_ST<ParticularKernel>::gather_Iss_() noexcept {

    const size_t t_ker_sz = this->up_kernel_->size();

    // Collect sub domains for each kernel
    this->Iss_tot_ = dcomplex(0.,0.);
    for (size_t i = 0; i < t_ker_sz; ++i) {
        dcomplex Iss_isum(0.0, 0.0);
        for (size_t kk = 0; kk < sub_figures_numb_k_(); ++kk) {
            Iss_isum += Isub_[i + t_ker_sz * kk];
        }
        this->Iss_[i] = Iss_isum * this->up_kernel_->precomputed_jacobian();
        this->Iss_tot_ += this->Iss_[i];
    }
}

//virtual
template <typename ParticularKernel>
void DirectfnAlgorithm_ST<ParticularKernel>::calc_n_functions_() noexcept {

    double Lam_b(0.0);
    setup_lamlim_(t_Eta_2_, Lam_b);  // get Lam_b

    const double half_lam_b = 0.5 * Lam_b;
    const double * const cz3 = this->up_z3_.get();
    const double * const cw3 = this->up_w3_.get();

    this->up_kerSummator_->nullify_Ilam_3();
    for (size_t n3_lam = 0; n3_lam < this->N3(); ++n3_lam) {

        const double t_Lam = half_lam_b * cz3[n3_lam] + half_lam_b;
        calc_a_functions_(t_Lam);
        this->up_kerSummator_->accumulate_Ilam_3(cw3[n3_lam]);
    }
    this->up_kerSummator_->multiply_Ilam_3(half_lam_b);
}

//virtual
template <typename ParticularKernel>
void DirectfnAlgorithm_ST<ParticularKernel>::calc_a_functions_(const double t_Lam) noexcept {

    // In the loop it takes Time(cos) x N4
    const double cos_PSI = cos_PSI_;
    const double sin_PSI = sin_PSI_;

    const double rho_b = t_Lam;
    const double half_rho_b = 0.5 * rho_b;

    const double etau_p_sub = t_Eta_2_;
    const double xiv_p_sub  = t_Lam * sin_PSI - tri_zero_qua_one_();  // 0 or 1

    double etau_p, xiv_p;
    pf_subtrg_crnt_(etau_p_sub, xiv_p_sub, etau_p, xiv_p); // in, in, out, out

    // Simplex coordinates for p
    double uvxi_p[3];
    st_make_simplex_(uvxi_p, etau_p, xiv_p);  // out[3], in, in
    this->up_kernel_->update_rp(uvxi_p);

    const double * const cz4 = this->up_z4_.get();
    const double * const cw4 = this->up_w4_.get();

    this->up_kerSummator_->nullify_Irho_4();
    for (size_t n4_rho = 0; n4_rho < this->N4(); ++n4_rho) {

        const double tmp_Rho = half_rho_b * cz4[n4_rho] + half_rho_b;

        const double etau_q_sub =  tmp_Rho * cos_PSI + etau_p_sub;
        const double xiv_q_sub  = -tmp_Rho * sin_PSI + xiv_p_sub;

        double etau_q, xiv_q;
        pf_subtrg_crnt_(etau_q_sub, xiv_q_sub, etau_q, xiv_q); // in, in, out, out

        // Simplex coordinates for q
        double uvxi_q[3];
        st_make_simplex_(uvxi_q, etau_q, xiv_q);  // out[3], in, in
        this->up_kernel_->update_rq(uvxi_q);

        // Prehashes the Greens function here.
        this->up_kernel_->precompute_rp_rq_data();
        this->up_kerSummator_->accumulate_Irho_4(cw4[n4_rho] * tmp_Rho);
    }
    this->up_kerSummator_->multiply_Irho_4(half_rho_b);
}

// Instantiation of the Triangular Constant Kernels
template class DirectfnAlgorithm_ST<TriangularKernel_Constant_ST>;
template class DirectfnAlgorithm_ST<TriangularKernel_Constant_EA>;
template class DirectfnAlgorithm_ST<TriangularKernel_Constant_VA>;
template class DirectfnAlgorithm_ST<TriangularKernel_RWG_WS>;

// Instantiation of the Quadrilateral Planar Kernels
template class DirectfnAlgorithm_ST<QuadrilateralKernel_PlanarScalar>;
template class DirectfnAlgorithm_ST<QuadrilateralKernel_PlanarVectorWS>;

template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class DirectfnAlgorithm_ST<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class DirectfnAlgorithm_ST<QuadrilateralKernel_CurvilinearScalar>;
template class DirectfnAlgorithm_ST<QuadrilateralKernel_CurvilinearVectorWS>;

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Triangular_ST<ParticularKernel>::Triangular_ST():
DirectfnAlgorithm_ST<ParticularKernel>(),
pf_subtrg_arr_{nullptr, nullptr, nullptr},
pf_psilim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_etalim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pmf_lamlim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pmf_lamlim_crnt_(nullptr) {

    // The kernel pointer as well as the kernel-array-summator pointer
    // has been allocated in the DirectfnInterface<ParticularKernel>

    // Here the only resources of the Triangular_ST are allocated.
    allocate_I_vars_();
    initialize_limit_fptrs_();
}

template <typename ParticularKernel>
Triangular_ST<ParticularKernel>::~Triangular_ST() {
}

//virtual
template <typename ParticularKernel>
string Triangular_ST<ParticularKernel>::name() const noexcept {
    return string("Triangular_ST");
}

template <typename ParticularKernel>
void Triangular_ST<ParticularKernel>::allocate_I_vars_() {

    const size_t t_ker_sz = this->up_kernel_->size();

    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_figures_3_() * t_ker_sz]);
    this->Isd_k_.reset(new dcomplex[sub_ranges_8_()* t_ker_sz]);
}

template <typename ParticularKernel>
void Triangular_ST<ParticularKernel>::initialize_limit_fptrs_() {

    // Array of templated function pointers to be used in loops when specialization is impossible
    // (at compile-time moment). That is why we need these wrappers.
    pf_subtrg_arr_[0] = triag_st_subtriangles<1>;
    pf_subtrg_arr_[1] = triag_st_subtriangles<2>;
    pf_subtrg_arr_[2] = triag_st_subtriangles<3>;

    // N1 psi limits initalization ptrs
    pf_psilim_arr_[0] = triag_st_psi_limits<1>;
    pf_psilim_arr_[1] = triag_st_psi_limits<2>;
    pf_psilim_arr_[2] = triag_st_psi_limits<3>;
    pf_psilim_arr_[3] = triag_st_psi_limits<4>;
    pf_psilim_arr_[4] = triag_st_psi_limits<5>;
    pf_psilim_arr_[5] = triag_st_psi_limits<6>;
    pf_psilim_arr_[6] = triag_st_psi_limits<7>;
    pf_psilim_arr_[7] = triag_st_psi_limits<8>;

    // N2 eta limits initialization ptrs
    pf_etalim_arr_[0] = triag_st_eta_limits<1>;
    pf_etalim_arr_[1] = triag_st_eta_limits<2>;
    pf_etalim_arr_[2] = triag_st_eta_limits<3>;
    pf_etalim_arr_[3] = triag_st_eta_limits<4>;
    pf_etalim_arr_[4] = triag_st_eta_limits<5>;
    pf_etalim_arr_[5] = triag_st_eta_limits<6>;
    pf_etalim_arr_[6] = triag_st_eta_limits<7>;
    pf_etalim_arr_[7] = triag_st_eta_limits<8>;

    // N3 lam limits initialization ptrs
    pmf_lamlim_arr_[0] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_1_;
    pmf_lamlim_arr_[1] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_2_;
    pmf_lamlim_arr_[2] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_3_;
    pmf_lamlim_arr_[3] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_4_;
    pmf_lamlim_arr_[4] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_5_;
    pmf_lamlim_arr_[5] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_6_;
    pmf_lamlim_arr_[6] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_7_;
    pmf_lamlim_arr_[7] = &Triangular_ST<ParticularKernel>::triag_st_lam_limits_8_;
}

//virtual
template <typename ParticularKernel>
size_t Triangular_ST<ParticularKernel>::sub_figures_numb_k_() const noexcept {
    return sub_figures_3_();
}

//virtual
template <typename ParticularKernel>
size_t Triangular_ST<ParticularKernel>::sub_ranges_numb_m_()  const noexcept {
    return sub_ranges_8_();
}

//virtual
template <typename ParticularKernel>
void Triangular_ST<ParticularKernel>::update_subtri_fptrs_(const size_t kk) noexcept {
    this->pf_subtrg_crnt_ = pf_subtrg_arr_[kk];
}

//virtual
template <typename ParticularKernel>
void Triangular_ST<ParticularKernel>::update_psilims_fptrs_N1_(const size_t m) noexcept {
    this->pf_psilim_crnt_ = pf_psilim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Triangular_ST<ParticularKernel>::update_etalims_fptrs_N2_(const size_t m) noexcept {
    this->pf_etalim_crnt_ = pf_etalim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Triangular_ST<ParticularKernel>::update_lamlims_fptrs_N3_(const size_t m) noexcept {
    this->pmf_lamlim_crnt_ = pmf_lamlim_arr_[m];
}

//virtual
template <typename ParticularKernel>
double Triangular_ST<ParticularKernel>::tri_zero_qua_one_() const noexcept {
    return 0.0;  // zero is proper value. It is not like
                 // "default" for nothing. For Triangular it is 0. Do not modify it.
}

//virtual
template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::st_make_simplex_(double chi[3], const double eta, const double xi) noexcept {
    triag_st_make_simplex(chi, eta, xi);  // out[3], in, in
}

//virtual
template <typename ParticularKernel>
void Triangular_ST<ParticularKernel>::setup_lamlim_(const double t_Eta_2, double & Lam_b) const noexcept {
    (this->*pmf_lamlim_crnt_)(t_Eta_2, Lam_b);
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_1_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = (double(1.0) - Eta) / this->cos_PSI_;
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_2_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = sqrt(double(3.0)) * (double(1.0) - Eta) / this->sin_PSI_;
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_3_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = sqrt(double(3.0)) * (double(1.0) - Eta) / this->sin_PSI_;
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_4_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = -(double(1.0) + Eta) / this->cos_PSI_;
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_5_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = (double(1.0) - Eta) / this->cos_PSI_;
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_6_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = sqrt(double(3.0)) * (double(1.0) + Eta) / this->sin_PSI_;
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_7_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = sqrt(double(3.0)) * (double(1.0) + Eta) / this->sin_PSI_;
}

template <typename ParticularKernel>
void  Triangular_ST<ParticularKernel>::triag_st_lam_limits_8_(const double Eta, double & Lam_b) const noexcept {
    Lam_b = -(double(1.0) + Eta) / this->cos_PSI_;
}


template class Triangular_ST<TriangularKernel_Constant_ST>;
template class Triangular_ST<TriangularKernel_Constant_EA>;
template class Triangular_ST<TriangularKernel_Constant_VA>;
template class Triangular_ST<TriangularKernel_RWG_WS>;

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Quadrilateral_ST<ParticularKernel>::Quadrilateral_ST():
DirectfnAlgorithm_ST<ParticularKernel>(),
cos_pi_ms_PSI_(0.0),
pf_subtrg_arr_{nullptr, nullptr, nullptr, nullptr},
pf_psilim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pf_etalim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pmf_lamlim_arr_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr},
pmf_lamlim_crnt_(nullptr) {

    // The kernel pointer as well as the kernel-array-summator pointer
    // has been allocated in the DirectfnInterface<ParticularKernel>

    // Here the memory for arrays is allocated according to kernel size.
    allocate_I_vars_();
    initialize_limit_fptrs_();
}

template <typename ParticularKernel>
Quadrilateral_ST<ParticularKernel>::~Quadrilateral_ST() {
}

//virtual
template <typename ParticularKernel>
string Quadrilateral_ST<ParticularKernel>::name() const noexcept {
    return string("Quadrilateral_ST");
}

//virtual
template <typename ParticularKernel>
void  Quadrilateral_ST<ParticularKernel>::precompute_psi_trigonomety_(const double t_Psi) noexcept {

    DirectfnAlgorithm_ST<ParticularKernel>::precompute_psi_trigonomety_(t_Psi);

    //tan_PSI_ms_05_pi_ = tan(Psi - 0.5 * M_PI);
    cos_pi_ms_PSI_ = cos(M_PI - t_Psi);
}

template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::allocate_I_vars_() {

    const size_t t_ker_sz = this->up_kernel_->size();

    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    allocate_I_vars_(t_ker_sz);
//    Iss_.reset(new dcomplex[t_ker_sz]);
//    Isub_.reset (new dcomplex[sub_figures_4_() * t_ker_sz]);
//    Isd_k_.reset(new dcomplex[sub_ranges_6_() * t_ker_sz]);
}

template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::allocate_I_vars_(const size_t t_ker_sz) noexcept {

    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_figures_4_() * t_ker_sz]);
    this->Isd_k_.reset(new dcomplex[sub_ranges_6_() * t_ker_sz]);
}

template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::initialize_limit_fptrs_() {

    // Array of templated function pointers to be used in loops when specialization is impossible
    // (at compile-time moment). That is why we need these wrappers.
    pf_subtrg_arr_[0] = quadr_st_subtriangles<1>;
    pf_subtrg_arr_[1] = quadr_st_subtriangles<2>;
    pf_subtrg_arr_[2] = quadr_st_subtriangles<3>;
    pf_subtrg_arr_[3] = quadr_st_subtriangles<4>;
    // N1
    pf_psilim_arr_[0] = quad_st_psi_limits<1>;
    pf_psilim_arr_[1] = quad_st_psi_limits<2>;
    pf_psilim_arr_[2] = quad_st_psi_limits<3>;
    pf_psilim_arr_[3] = quad_st_psi_limits<4>;
    pf_psilim_arr_[4] = quad_st_psi_limits<5>;
    pf_psilim_arr_[5] = quad_st_psi_limits<6>;
    // N2
    pf_etalim_arr_[0] = quad_st_eta_limits<16>;
    pf_etalim_arr_[1] = quad_st_eta_limits<2>;
    pf_etalim_arr_[2] = quad_st_eta_limits<3>;
    pf_etalim_arr_[3] = quad_st_eta_limits<4>;
    pf_etalim_arr_[4] = quad_st_eta_limits<5>;
    pf_etalim_arr_[5] = quad_st_eta_limits<16>;
    // N3
    pmf_lamlim_arr_[0] = &Quadrilateral_ST::quad_st_lamlims_12_;
    pmf_lamlim_arr_[1] = &Quadrilateral_ST::quad_st_lamlims_12_;
    pmf_lamlim_arr_[2] = &Quadrilateral_ST::quad_st_lamlims_34_;
    pmf_lamlim_arr_[3] = &Quadrilateral_ST::quad_st_lamlims_34_;
    pmf_lamlim_arr_[4] = &Quadrilateral_ST::quad_st_lamlims_56_;
    pmf_lamlim_arr_[5] = &Quadrilateral_ST::quad_st_lamlims_56_;
}

//virtual
template <typename ParticularKernel>
size_t Quadrilateral_ST<ParticularKernel>::sub_figures_numb_k_() const noexcept {
    return sub_figures_4_();
}

//virtual
template <typename ParticularKernel>
size_t Quadrilateral_ST<ParticularKernel>::sub_ranges_numb_m_()  const noexcept {
    return sub_ranges_6_();
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::update_subtri_fptrs_(const size_t kk) noexcept {
    this->pf_subtrg_crnt_ = pf_subtrg_arr_[kk];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::update_psilims_fptrs_N1_(const size_t m) noexcept {
    this->pf_psilim_crnt_ = pf_psilim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::update_etalims_fptrs_N2_(const size_t m) noexcept {
    this->pf_etalim_crnt_ = pf_etalim_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::update_lamlims_fptrs_N3_(const size_t m) noexcept {
    this->pmf_lamlim_crnt_ = pmf_lamlim_arr_[m];
}

//virtual
template <typename ParticularKernel>
double Quadrilateral_ST<ParticularKernel>::tri_zero_qua_one_() const noexcept {
    return 1.0;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::st_make_simplex_(double chi[3], const double eta, const double xi) noexcept {

    chi[0] = eta;
    chi[1] = xi;
    chi[2] = 0.0;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::setup_lamlim_(const double t_Eta_2, double & Lam_b) const noexcept {
    (this->*pmf_lamlim_crnt_)(t_Eta_2, Lam_b);
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::quad_st_lamlims_12_(const double U_inp, double & Lam_b) const noexcept {
    Lam_b = (1.0 - U_inp) / this->cos_PSI_;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::quad_st_lamlims_34_(const double , double & Lam_b) const noexcept {
    Lam_b = 2.0 / this->sin_PSI_;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST<ParticularKernel>::quad_st_lamlims_56_(const double U_inp, double & Lam_b) const noexcept {
    Lam_b = (1.0 + U_inp) / cos_pi_ms_PSI_;
}


// Instantiation of the Quadrilateral Constant Kernels
template class Quadrilateral_ST<QuadrilateralKernel_PlanarScalar>;
template class Quadrilateral_ST<QuadrilateralKernel_PlanarVectorWS>;
template class Quadrilateral_ST<QuadrilateralKernel_PlanarVectorSS>;

template class Quadrilateral_ST<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class Quadrilateral_ST<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class Quadrilateral_ST<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class Quadrilateral_ST<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class Quadrilateral_ST<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class Quadrilateral_ST<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class Quadrilateral_ST<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class Quadrilateral_ST<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class Quadrilateral_ST<QuadrilateralKernel_CurvilinearScalar>;
template class Quadrilateral_ST<QuadrilateralKernel_CurvilinearVectorWS>;
template class Quadrilateral_ST<QuadrilateralKernel_CurvilinearVectorSS>;

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file




