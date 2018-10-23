#include <iostream>
#include <string>
#include "directfn_algorithm_va.h"
#include "directfn_kernel_tri.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"
#include "directfn_kernel_quad_voxhenry.h"

using  std::cout;
using  std::endl;
using  std::string;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template  <typename ParticularKernel>
DirectfnAlgorithm_VA<ParticularKernel>::DirectfnAlgorithm_VA():
DirectfnInterface<ParticularKernel>(),
pf_theta_p_lim_1_crnt_(nullptr),
pf_theta_q_lim_2_crnt_(nullptr),
pf_get_Lp_1_crnt_(nullptr),
pf_get_Lq_2_crnt_(nullptr),
pf_crnt_lam_maxlim_(nullptr),
Isub_(nullptr),
sin_Theta_p_1_(0.0),
cos_Theta_p_1_(0.0),
sin_Theta_q_2_(0.0),
cos_Theta_q_2_(0.0),
psi_a_(0.0),
psi_b_(0.0),
sin_Psi_3_(0.0),
cos_Psi_3_(0.0),
Lp_1_(0.0),
Lq_2_(0.0) {
}

template  <typename ParticularKernel>
DirectfnAlgorithm_VA<ParticularKernel>::~DirectfnAlgorithm_VA() {
}

//virtual
template  <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::do_I_surface_surface_() {

    // Go throw subranges
    for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {

        // Update limits f-pointers for a given m'th subrange
        update_theta_p_lims_N1_fptr_(m);
        update_theta_q_lims_N2_fptr_(m);
        update_Lp1_fptr_(m);
        update_Lq2_fptr_(m);

        double theta_p_A = 0.0, theta_p_B = 0.0;
        pf_theta_p_lim_1_crnt_(theta_p_A, theta_p_B);
        const double half_theta_p_a_plus_b = 0.5 * (theta_p_A + theta_p_B);
        const double half_theta_p_b_mnus_a = 0.5 * (theta_p_B - theta_p_A);

        const double * const cz1 = this->up_z1_.get();
        const double * const cw1 = this->up_w1_.get();

        this->up_kerSummator_->nullify_Ipsi_1();
        for (size_t n_theta_p_1 = 0; n_theta_p_1 < this->N1(); ++n_theta_p_1) {

            const double Theta_p1 = half_theta_p_b_mnus_a * cz1[n_theta_p_1] + half_theta_p_a_plus_b;
            update_theta_p1_trigonometry_(Theta_p1);
            Lp_1_ = pf_get_Lp_1_crnt_(cos_Theta_p_1_, sin_Theta_p_1_);

            double theta_q_A = 0.0, theta_q_B = 0.0;
            pf_theta_q_lim_2_crnt_(theta_q_A, theta_q_B);
            const double half_theta_q_a_plus_b = 0.5 * (theta_q_A + theta_q_B);
            const double half_theta_q_b_mnus_a = 0.5 * (theta_q_B - theta_q_A);

            const double * const cz2 = this->up_z2_.get();
            const double * const cw2 = this->up_w2_.get();

            this->up_kerSummator_->nullify_Ieta_2();
            for (size_t n_theta_q_2 = 0; n_theta_q_2 < this->N2(); ++n_theta_q_2) {

                const double Theta_q2 = half_theta_q_b_mnus_a * cz2[n_theta_q_2] + half_theta_q_a_plus_b;
                update_theta_q2_trigonometry_(Theta_q2);
                Lq_2_ = pf_get_Lq_2_crnt_(cos_Theta_q_2_, sin_Theta_q_2_);
                // Lp1 and Lq2 are ready now:
                const double atan_Lq_to_Lp = atan(Lq_2_ / Lp_1_);
                // 1
                define_psi_a_psi_b_lamlim3_case_I_1_(atan_Lq_to_Lp);
                calc_I_psi_N3_();
                this->up_kerSummator_->accumulate_Ieta_2(cw2[n_theta_q_2]);
                // 2
                define_psi_a_psi_b_lamlim3_case_II_2_(atan_Lq_to_Lp);
                calc_I_psi_N3_();
                this->up_kerSummator_->accumulate_Ieta_2(cw2[n_theta_q_2]);
            } // N2 for
            this->up_kerSummator_->multiply_Ieta_2(half_theta_q_b_mnus_a);
            this->up_kerSummator_->accumulate_Ipsi_1(cw1[n_theta_p_1]);
        } // N1 for
        this->up_kerSummator_->multiply_Ipsi_1(half_theta_p_b_mnus_a);
        this->up_kerSummator_->assign_Ipsi_1_to( &Isub_[this->up_kernel_->size() * m] );
    } // 1 or 4  ==  m for triangles or quadrilateral
    gather_Iss_();
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::gather_Iss_() noexcept {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Collect sub ranges for each kernel
    this->Iss_tot_ = dcomplex(0.,0.);
    for (size_t i = 0; i < t_ker_sz; ++i) {
        dcomplex Iss_isum(0.0, 0.0);
        for (size_t m = 0; m < sub_ranges_numb_m_(); ++m) {
            Iss_isum += this->Isub_[i + t_ker_sz * m];
        }
        this->Iss_[i] = Iss_isum * this->up_kernel_->precomputed_jacobian();
        this->Iss_tot_ += this->Iss_[i];
    }
}

//virtual
template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::update_theta_p1_trigonometry_(const double Theta_p1) noexcept {

    sin_Theta_p_1_ = sin(Theta_p1);
    cos_Theta_p_1_ = cos(Theta_p1);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::update_theta_q2_trigonometry_(const double Theta_q2) noexcept {

    sin_Theta_q_2_ = sin(Theta_q2);
    cos_Theta_q_2_ = cos(Theta_q2);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::define_psi_a_psi_b_lamlim3_case_I_1_ (const double atan_Lq_to_Lp) noexcept {

    psi_a_ = 0.0;
    psi_b_ = atan_Lq_to_Lp;
    pf_crnt_lam_maxlim_ = &DirectfnAlgorithm_VA<ParticularKernel>::lam_lim_b_1_;
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::define_psi_a_psi_b_lamlim3_case_II_2_(const double atan_Lq_to_Lp) noexcept {

    psi_a_ = atan_Lq_to_Lp;
    psi_b_ = M_PI / 2.0;
    pf_crnt_lam_maxlim_ = &DirectfnAlgorithm_VA<ParticularKernel>::lam_lim_b_2_;
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::calc_I_psi_N3_() noexcept {

    const double psi_plus_ab = 0.5 * (psi_b_ + psi_a_);
    const double psi_mnus_ab = 0.5 * (psi_b_ - psi_a_);

    const double * const cz3 = this->up_z3_.get();
    const double * const cw3 = this->up_w3_.get();

    this->up_kerSummator_->nullify_Ilam_3();
    for (size_t n3_psi = 0; n3_psi < this->N3(); ++n3_psi) {

        const double t_Psi = psi_mnus_ab * cz3[n3_psi] + psi_plus_ab;
        update_psi_3_trigonometry_(t_Psi);

        calc_I_rho_4_();
        this->up_kerSummator_->accumulate_Ilam_3(cw3[n3_psi]);
    }
    this->up_kerSummator_->multiply_Ilam_3(psi_mnus_ab);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::update_psi_3_trigonometry_(const double t_Psi_3) noexcept {

    cos_Psi_3_ = cos(t_Psi_3);
    sin_Psi_3_ = sin(t_Psi_3);
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::calc_I_rho_4_() noexcept {

    const double J_lam = 0.5 * (this->*pf_crnt_lam_maxlim_)();

    const double * const cz4 = this->up_z4_.get();
    const double * const cw4 = this->up_w4_.get();

    this->up_kerSummator_->nullify_Irho_4();
    for (size_t n4_lam = 0; n4_lam < this->N4(); ++n4_lam) {

        const double t_Lam = J_lam * (cz4[n4_lam] + 1.0);
        double uvxi_p[3], uvxi_q[3];
        calc_pq_simplex_(uvxi_p, uvxi_q, t_Lam); // out, out, in

        this->up_kernel_->update_rp(uvxi_p);
        this->up_kernel_->update_rq(uvxi_q);
        // Precaches the Greens function here.
        this->up_kernel_->precompute_rp_rq_data();
        this->up_kerSummator_->accumulate_Irho_4(cw4[n4_lam] * t_Lam * t_Lam * t_Lam);
    }
    this->up_kerSummator_->multiply_Irho_4(sin_Psi_3_ * cos_Psi_3_ * J_lam);
}

template <typename ParticularKernel>
double DirectfnAlgorithm_VA<ParticularKernel>::lam_lim_b_1_() const noexcept {
    return Lp_1_ / cos_Psi_3_;
}

template <typename ParticularKernel>
double DirectfnAlgorithm_VA<ParticularKernel>::lam_lim_b_2_() const noexcept {
    return Lq_2_ / sin_Psi_3_;
}

template <typename ParticularKernel>
void DirectfnAlgorithm_VA<ParticularKernel>::calc_pq_simplex_(double xi_p_out[3], double xi_q_out[3],
                                            const double t_Lam) noexcept {

    const double u_eta_p = t_Lam * cos_Psi_3_ * cos_Theta_p_1_ - 1.0;
    const double v_xi_p  = t_Lam * cos_Psi_3_ * sin_Theta_p_1_ - unity4_zero3_();

    va_make_simplex(xi_p_out, u_eta_p, v_xi_p);

    const double u_eta_q = t_Lam * sin_Psi_3_ * cos_Theta_q_2_ - 1.0;
    const double v_xi_q  = t_Lam * sin_Psi_3_ * sin_Theta_q_2_ - unity4_zero3_();

    va_make_simplex(xi_q_out, u_eta_q, v_xi_q);
}

// Instantiation of the Triangular Constant Kernels
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_ST>;
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_EA>;
template class DirectfnAlgorithm_VA<TriangularKernel_Constant_VA>;

template class DirectfnAlgorithm_VA<TriangularKernel_RWG_WS>;
template class DirectfnAlgorithm_VA<TriangularKernel_RWG_SS>;
template class DirectfnAlgorithm_VA<TriangularKernel_nxRWG_SS>;

// Instantiation of the Quadrilateral Planar Kernels
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarScalar>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarVectorWS>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_PlanarVectorSS>;
// .. nxRWG?

template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class DirectfnAlgorithm_VA<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearScalar>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearVectorWS>;
template class DirectfnAlgorithm_VA<QuadrilateralKernel_CurvilinearVectorSS>;

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Triangular_VA<ParticularKernel>::Triangular_VA():
DirectfnAlgorithm_VA<ParticularKernel>() {

    // The kernel pointer as well as the kernel-array-summator pointer
    // has been allocated in the DirectfnInterface<ParticularKernel>

    // Here the memory for arrays is allocated according to kernel size.
    allocate_I_vars_();
    initialize_limit_fptrs_();
}

template <typename ParticularKernel>
Triangular_VA<ParticularKernel>::~Triangular_VA() {
}

//virtual
template <typename ParticularKernel>
string Triangular_VA<ParticularKernel>::name() const noexcept {
    return string("Triangular_VA_Constant");
}

template <typename ParticularKernel>
void Triangular_VA<ParticularKernel>::allocate_I_vars_() {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_figures_1_() * t_ker_sz]);
}

template <typename ParticularKernel>
void Triangular_VA<ParticularKernel>::initialize_limit_fptrs_() noexcept {

    this->pf_theta_p_lim_1_crnt_ = triag_va_theta_p_limits;
    this->pf_theta_q_lim_2_crnt_ = triag_va_theta_q_limits;

    this->pf_get_Lp_1_crnt_ = triag_va_Lp1;
    this->pf_get_Lq_2_crnt_ = triag_va_Lq2;
}

//virtual
template <typename ParticularKernel>
size_t Triangular_VA<ParticularKernel>::sub_ranges_numb_m_()  const noexcept {
    return sub_figures_1_();
}

//virtual
template <typename ParticularKernel>
void Triangular_VA<ParticularKernel>::update_theta_p_lims_N1_fptr_(const size_t ) noexcept {
    // function pointer  pf_theta_p_lim_1_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel>
void Triangular_VA<ParticularKernel>::update_theta_q_lims_N2_fptr_(const size_t ) noexcept {
    // function pointer pf_theta_q_lim_2_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel>
void Triangular_VA<ParticularKernel>::update_Lp1_fptr_(const size_t ) noexcept {
    // function pointer pf_get_Lp_1_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel>
void Triangular_VA<ParticularKernel>::update_Lq2_fptr_(const size_t ) noexcept {
    // function pointer pf_get_Lq_2_crnt_
    // have been already setup in the constructor. No need to reset.
}

//virtual
template <typename ParticularKernel>
double Triangular_VA<ParticularKernel>::unity4_zero3_() const noexcept {
    return 0.0;
}

//virtual
template <typename ParticularKernel>
void Triangular_VA<ParticularKernel>::va_make_simplex(double uvxi_p_out[3], const double u_eta_p,
                                                  const double v_xi_p) const noexcept {
    triag_st_make_simplex(uvxi_p_out, u_eta_p, v_xi_p);   // out[3], in, in
}


template class Triangular_VA<TriangularKernel_Constant_ST>;
template class Triangular_VA<TriangularKernel_Constant_EA>;
template class Triangular_VA<TriangularKernel_Constant_VA>;

template class Triangular_VA<TriangularKernel_RWG_WS>;
template class Triangular_VA<TriangularKernel_RWG_SS>;
template class Triangular_VA<TriangularKernel_nxRWG_SS>;

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Quadrilateral_VA<ParticularKernel>::Quadrilateral_VA():
DirectfnInterface<ParticularKernel>(),
pF_theta_p_lim1_arr_{nullptr, nullptr, nullptr, nullptr},
pF_theta_q_lim2_arr_{nullptr, nullptr, nullptr, nullptr},
pf_Lp_1_arr_{nullptr, nullptr, nullptr, nullptr},
pf_Lq_2_arr_{nullptr, nullptr, nullptr, nullptr} {

    // Initializtion of particular kernels is done in the
    // DirectfnInterface<ParticularKernel>
    allocate_I_vars_();
    initialize_limit_fptrs_();
}

template <typename ParticularKernel>
Quadrilateral_VA<ParticularKernel>::~Quadrilateral_VA() {
}

//virtual
template <typename ParticularKernel>
string Quadrilateral_VA<ParticularKernel>::name() const noexcept {
    return string("Quadrilateral_VA");
}

template <typename ParticularKernel>
void Quadrilateral_VA<ParticularKernel>::allocate_I_vars_() {

    const size_t t_ker_sz = this->up_kernel_->size();
    // Allocates the memory according to the bt-type (the kernel length may vary
    // depending on type of Basis-Testing functions)
    // The kernel vector is garanteed to have finite length
    this->Iss_.reset(new dcomplex[t_ker_sz]);
    this->Isub_.reset (new dcomplex[sub_figures_4_() * t_ker_sz]);
}

template <typename ParticularKernel>
void Quadrilateral_VA<ParticularKernel>::initialize_limit_fptrs_() noexcept {
    // N1
    pF_theta_p_lim1_arr_[0] = quad_va_theta_p_limits<0>;
    pF_theta_p_lim1_arr_[1] = quad_va_theta_p_limits<1>;
    pF_theta_p_lim1_arr_[2] = quad_va_theta_p_limits<2>;
    pF_theta_p_lim1_arr_[3] = quad_va_theta_p_limits<3>;
    // N2
    pF_theta_q_lim2_arr_[0] = quad_va_theta_q_limits<0>;
    pF_theta_q_lim2_arr_[1] = quad_va_theta_q_limits<1>;
    pF_theta_q_lim2_arr_[2] = quad_va_theta_q_limits<2>;
    pF_theta_q_lim2_arr_[3] = quad_va_theta_q_limits<3>;
    // N3
    pf_Lp_1_arr_[0] = quad_va_Lp1<0>;
    pf_Lp_1_arr_[1] = quad_va_Lp1<1>;
    pf_Lp_1_arr_[2] = quad_va_Lp1<2>;
    pf_Lp_1_arr_[3] = quad_va_Lp1<3>;
    // N4
    pf_Lq_2_arr_[0] = quad_va_Lq2<0>;
    pf_Lq_2_arr_[1] = quad_va_Lq2<1>;
    pf_Lq_2_arr_[2] = quad_va_Lq2<2>;
    pf_Lq_2_arr_[3] = quad_va_Lq2<3>;
}

//virtual
template <typename ParticularKernel>
size_t Quadrilateral_VA<ParticularKernel>::sub_ranges_numb_m_() const noexcept {
    return sub_figures_4_();
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_VA<ParticularKernel>::update_theta_p_lims_N1_fptr_(const size_t  m) noexcept {
    this->pf_theta_p_lim_1_crnt_ = pF_theta_p_lim1_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_VA<ParticularKernel>::update_theta_q_lims_N2_fptr_(const size_t  m) noexcept {
    this->pf_theta_q_lim_2_crnt_ = pF_theta_q_lim2_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_VA<ParticularKernel>::update_Lp1_fptr_(const size_t m) noexcept {
    this->pf_get_Lp_1_crnt_ = pf_Lp_1_arr_[m];
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_VA<ParticularKernel>::update_Lq2_fptr_(const size_t m) noexcept {
    this->pf_get_Lq_2_crnt_ = pf_Lq_2_arr_[m];
}

//virtual
template <typename ParticularKernel>
double Quadrilateral_VA<ParticularKernel>::unity4_zero3_() const noexcept {
    return 1.0;
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_VA<ParticularKernel>::va_make_simplex(double uvxi_p_out[3], const double u_eta_p,
                                                               const double v_xi_p) const noexcept {
    uvxi_p_out[0] = u_eta_p;
    uvxi_p_out[1] = v_xi_p;
    uvxi_p_out[2] = 0.0;
}


// Instantiation of the Quadrilateral Constant Kernels
template class Quadrilateral_VA<QuadrilateralKernel_PlanarScalar>;
template class Quadrilateral_VA<QuadrilateralKernel_PlanarVectorWS>;
template class Quadrilateral_VA<QuadrilateralKernel_PlanarVectorSS>;

template class Quadrilateral_VA<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class Quadrilateral_VA<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class Quadrilateral_VA<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class Quadrilateral_VA<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class Quadrilateral_VA<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class Quadrilateral_VA<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class Quadrilateral_VA<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class Quadrilateral_VA<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearScalar>;
template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorWS>;
template class Quadrilateral_VA<QuadrilateralKernel_CurvilinearVectorSS>;

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file

