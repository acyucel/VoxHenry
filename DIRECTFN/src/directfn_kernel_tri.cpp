#include <iostream>
#include <iomanip>
#include "directfn_kernel_tri.h"
#include "directfn_interface.h"
#include "directfn_greenfunc.h"

using  std::cout;
using  std::cerr;
using  std::endl;
using  std::setprecision;
using  std::setw;
using  std::to_string;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

TriangularKernel::TriangularKernel():
AbstractKernel(),
rp1_{0.0, 0.0, 0.0},
rp2_{0.0, 0.0, 0.0},
rp3_{0.0, 0.0, 0.0},
rq1_{0.0, 0.0, 0.0},
rq2_{0.0, 0.0, 0.0},
rq3_{0.0, 0.0, 0.0} {

}

//virtual
TriangularKernel::~TriangularKernel() {
}

//virtual
void TriangularKernel::update_rp(const double uvxi_p[3]) noexcept {
    for (int i = 0; i < 3; ++i) {
        rp_crnt_[i] = uvxi_p[0] * rp1_[i] + uvxi_p[1] * rp2_[i] + uvxi_p[2] * rp3_[i];
    }
}

//virtual
void TriangularKernel::update_rq(const double uvxi_q[3]) noexcept {
    for (int i = 0; i < 3; ++i) {
        rq_crnt_[i] = uvxi_q[0] * rq1_[i] + uvxi_q[1] * rq2_[i] + uvxi_q[2] * rq3_[i];
    }
}

//virtual
void TriangularKernel::set(const SingularContour3xn &  contour_xpts) noexcept {

    if (3 == contour_xpts.length()) {
        set_3_pts_(contour_xpts);
        return;
    }
    if (4 == contour_xpts.length()) {
        set_4_pts_(contour_xpts);
        return;
    }
    if (5 == contour_xpts.length()) {
        set_5_pts_(contour_xpts);
        return;
    }
    cerr << "Too few/many vertexes to setup quadrilateral element" << endl;
}

//virtual
void  TriangularKernel::set_3_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r1_in = contour_xpts(0);
    const double * r2_in = contour_xpts(1);
    const double * r3_in = contour_xpts(2);

    for (size_t i = 0; i < 3; ++i) {

        this->rp1_[i] = r1_in[i];
        this->rp2_[i] = r2_in[i];
        this->rp3_[i] = r3_in[i];

        this->rq1_[i] = r1_in[i];
        this->rq2_[i] = r2_in[i];
        this->rq3_[i] = r3_in[i];
    }
    // lp[3], lq[3], np_[3] if needed
    precompute_vertex_dependent_data_();
}

//virtual
void  TriangularKernel::set_4_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r1 = contour_xpts(0);
    const double * r2 = contour_xpts(1);
    const double * r3 = contour_xpts(2);
    const double * r4 = contour_xpts(3);

    for (size_t i = 0; i < 3; ++i) {

        this->rp1_[i] = r1[i];
        this->rp2_[i] = r2[i];
        this->rp3_[i] = r3[i];

        this->rq1_[i] = r2[i];
        this->rq2_[i] = r1[i];
        this->rq3_[i] = r4[i];
    }
    // lp[3], lq[3], np_[3] if needed
    precompute_vertex_dependent_data_();
}

//virtual
void  TriangularKernel::set_5_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r1 = contour_xpts(0);
    const double * r2 = contour_xpts(1);
    const double * r3 = contour_xpts(2);
    const double * r4 = contour_xpts(3);
    const double * r5 = contour_xpts(4);

    for (size_t i = 0; i < 3; ++i) {

        this->rp1_[i] = r1[i];
        this->rp2_[i] = r2[i];
        this->rp3_[i] = r3[i];

        this->rq1_[i] = r1[i];
        this->rq2_[i] = r4[i];
        this->rq3_[i] = r5[i];
    }
    // lp[3], lq[3], np_[3] if needed
    precompute_vertex_dependent_data_();
}

//virtual
void TriangularKernel::debug_print() const noexcept {
    // TO DO:
    cout << " Warning! Empty  void TriangularKernel::debug_print() const noexcept " << endl;
}

//virtual
void TriangularKernel::precompute_rp_derivatives_() noexcept {
    // Empty definition for non abstract class.
    // Will be implemented later with curvilinear elements.
}

//virtual
void TriangularKernel::precompute_rq_derivatives_() noexcept {
    // Empty definition for non abstract class.
    // Will be implemented later with curvilinear elements.
}

//virtual
double TriangularKernel::calc_local_Jacobian_() noexcept {
    // Empty definition for non abstract class.
    // Will be implemented later with curvilinear elements.
    return 1.0;
}

///////////////////////////////////////////////////////////////////////////////

TriangularKernel_Constant::TriangularKernel_Constant():
TriangularKernel(),
jacobian_(0.0) {
    up_green_func_.reset(new HelmgolzGreenFunc());
}

//virtual
TriangularKernel_Constant::~TriangularKernel_Constant() {
}

//virtual
size_t TriangularKernel_Constant::size() const noexcept {
    return TriangularKernel_Constant::constexpr_size();
}

//virtual
double TriangularKernel_Constant::precomputed_jacobian() const noexcept {
    return jacobian_;
}

//virtual
void TriangularKernel_Constant::precompute_rp_rq_dependent_data_() noexcept {
}

//virtual
dcomplex TriangularKernel_Constant::specific_value_(const size_t ) const noexcept {
    return up_green_func_->value();
}

///----------------------------------------------------------------------------
// virtual
void TriangularKernel_Constant_ST::precompute_vertex_dependent_data_() noexcept {

    double l_2[3], l_3[3];
    for (size_t i = 0; i < 3; ++i) {
        l_2[i] = rp3_[i] - rp1_[i];
        l_3[i] = rp1_[i] - rp2_[i];
    }
    double Ap_v[3];
    vector_cross(l_3, l_2, Ap_v);
    jacobian_ = vector_dot(Ap_v, Ap_v) / double(12.0);
}

///----------------------------------------------------------------------------
//virtual
void TriangularKernel_Constant_EA::precompute_vertex_dependent_data_() noexcept {

    double lp_2[3], lp_3[3], lq_2[3], lq_3[3];
    for (size_t i = 0; i < 3; ++i) {
        lp_2[i] = rp3_[i] - rp1_[i];
        lp_3[i] = rp2_[i] - rp1_[i];
        lq_2[i] = rq3_[i] - rq1_[i];
        lq_3[i] = rq2_[i] - rq1_[i];
    }

    double Ap_v[3], Aq_v[3];
    vector_cross(lp_3, lp_2, Ap_v);
    vector_cross(lq_3, lq_2, Aq_v);

    jacobian_ = sqrt(vector_dot(Ap_v, Ap_v)) * sqrt(vector_dot(Aq_v, Aq_v)) / 12.0;
}

///----------------------------------------------------------------------------
//virtual
void TriangularKernel_Constant_VA::precompute_vertex_dependent_data_() noexcept {

    double lp_2[3], lp_3[3], lq_2[3], lq_3[3];
    for (size_t i = 0; i < 3; ++i) {
        lp_2[i] = rp3_[i] - rp1_[i];
        lp_3[i] = rp2_[i] - rp1_[i];
        lq_2[i] = rq3_[i] - rq1_[i];
        lq_3[i] = rq2_[i] - rq1_[i];
    }

    double Ap_v[3], Aq_v[3];
    vector_cross(lp_3, lp_2, Ap_v);
    vector_cross(lq_3, lq_2, Aq_v);

    jacobian_ = sqrt(vector_dot(Ap_v, Ap_v)) * sqrt(vector_dot(Aq_v, Aq_v)) / 12.0;
}

///////////////////////////////////////////////////////////////////////////////

/*! Specialization must be defined before call according to C++ Standard */
template <> dcomplex TriangularKernel_RWG::value_ind_<0>() const noexcept {

    cp_f_i_ = f_1_;
    cp_g_j_ = g_1_;

    lp_i_ = lp_[0];
    lq_j_ = lq_[0];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<1>() const noexcept {

    cp_f_i_ = f_1_;
    cp_g_j_ = g_2_;

    lp_i_ = lp_[0];
    lq_j_ = lq_[1];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<2>() const noexcept {

    cp_f_i_ = f_1_;
    cp_g_j_ = g_3_;

    lp_i_ = lp_[0];
    lq_j_ = lq_[2];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<3>() const noexcept {

    cp_f_i_ = f_2_;
    cp_g_j_ = g_1_;

    lp_i_ = lp_[1];
    lq_j_ = lq_[0];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<4>() const noexcept {

    cp_f_i_ = f_2_;
    cp_g_j_ = g_2_;

    lp_i_ = lp_[1];
    lq_j_ = lq_[1];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<5>() const noexcept {

    cp_f_i_ = f_2_;
    cp_g_j_ = g_3_;

    lp_i_ = lp_[1];
    lq_j_ = lq_[2];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<6>() const noexcept {

    cp_f_i_ = f_3_;
    cp_g_j_ = g_1_;

    lp_i_ = lp_[2];
    lq_j_ = lq_[0];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<7>() const noexcept {

    cp_f_i_ = f_3_;
    cp_g_j_ = g_2_;

    lp_i_ = lp_[2];
    lq_j_ = lq_[1];

    return this->rwg_value_();
}

template <> dcomplex TriangularKernel_RWG::value_ind_<8>() const noexcept {

    cp_f_i_ = f_3_;
    cp_g_j_ = g_3_;

    lp_i_ = lp_[2];
    lq_j_ = lq_[2];

    return this->rwg_value_();
}

TriangularKernel_RWG::TriangularKernel_RWG():
TriangularKernel(),
lp_{0.0, 0.0, 0.0},
lq_{0.0, 0.0, 0.0},
f_1_{0.0, 0.0, 0.0},
f_2_{0.0, 0.0, 0.0},
f_3_{0.0, 0.0, 0.0},
g_1_{0.0, 0.0, 0.0},
g_2_{0.0, 0.0, 0.0},
g_3_{0.0, 0.0, 0.0},
cp_f_i_(nullptr),
cp_g_j_(nullptr),
lp_i_(0.0),
lq_j_(0.0),
pmf_rwg_val_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr} {

    pmf_rwg_val_[0] = &TriangularKernel_RWG::value_ind_<0>;
    pmf_rwg_val_[1] = &TriangularKernel_RWG::value_ind_<1>;
    pmf_rwg_val_[2] = &TriangularKernel_RWG::value_ind_<2>;
    pmf_rwg_val_[3] = &TriangularKernel_RWG::value_ind_<3>;
    pmf_rwg_val_[4] = &TriangularKernel_RWG::value_ind_<4>;
    pmf_rwg_val_[5] = &TriangularKernel_RWG::value_ind_<5>;
    pmf_rwg_val_[6] = &TriangularKernel_RWG::value_ind_<6>;
    pmf_rwg_val_[7] = &TriangularKernel_RWG::value_ind_<7>;
    pmf_rwg_val_[8] = &TriangularKernel_RWG::value_ind_<8>;
}

TriangularKernel_RWG::~TriangularKernel_RWG() {
}

size_t TriangularKernel_RWG::size() const noexcept {
    return TriangularKernel_RWG::constexpr_size();
}

double TriangularKernel_RWG::precomputed_jacobian() const noexcept {
    return 1.0;
}

//virtual
void TriangularKernel_RWG::precompute_vertex_dependent_data_() noexcept {

    const double  ev_1p[3] = {rp3_[0] - rp2_[0], rp3_[1] - rp2_[1], rp3_[2] - rp2_[2]};
    const double  ev_2p[3] = {rp1_[0] - rp3_[0], rp1_[1] - rp3_[1], rp1_[2] - rp3_[2]};
    const double  ev_3p[3] = {rp2_[0] - rp1_[0], rp2_[1] - rp1_[1], rp2_[2] - rp1_[2]};

    const double  ev_1q[3] = {rq3_[0] - rq2_[0], rq3_[1] - rq2_[1], rq3_[2] - rq2_[2]};
    const double  ev_2q[3] = {rq1_[0] - rq3_[0], rq1_[1] - rq3_[1], rq1_[2] - rq3_[2]};
    const double  ev_3q[3] = {rq2_[0] - rq1_[0], rq2_[1] - rq1_[1], rq2_[2] - rq1_[2]};

    lp_[0] = sqrt(vector_dot(ev_1p, ev_1p));
    lp_[1] = sqrt(vector_dot(ev_2p, ev_2p));
    lp_[2] = sqrt(vector_dot(ev_3p, ev_3p));

    lq_[0] = sqrt(vector_dot(ev_1q, ev_1q));
    lq_[1] = sqrt(vector_dot(ev_2q, ev_2q));
    lq_[2] = sqrt(vector_dot(ev_3q, ev_3q));
}

//virtual
void TriangularKernel_RWG::precompute_rp_rq_dependent_data_()  noexcept {

    for (size_t i = 0; i < 3; ++i) {

        f_1_[i] = rp_crnt_[i] - rp1_[i];
        f_2_[i] = rp_crnt_[i] - rp2_[i];
        f_3_[i] = rp_crnt_[i] - rp3_[i];

        g_1_[i] = rq_crnt_[i] - rq1_[i];
        g_2_[i] = rq_crnt_[i] - rq2_[i];
        g_3_[i] = rq_crnt_[i] - rq3_[i];
    }
}

//virtual
dcomplex TriangularKernel_RWG::specific_value_(const size_t index) const noexcept {
    return (this->*pmf_rwg_val_[index])() * up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////

TriangularKernel_RWG_WS::TriangularKernel_RWG_WS():
TriangularKernel_RWG() {
    up_green_func_.reset(new HelmgolzGreenFunc());
}

TriangularKernel_RWG_WS::~TriangularKernel_RWG_WS() {
}

//virtual
dcomplex TriangularKernel_RWG_WS::rwg_value_() const noexcept {
    return  vector_dot(cp_f_i_, cp_g_j_) * lp_i_ * lq_j_ / 12.0;
}

///////////////////////////////////////////////////////////////////////////////

TriangularKernel_RWG_SS::TriangularKernel_RWG_SS():
TriangularKernel_RWG() {
    up_green_func_.reset(new GradHelmgolzGreenFunc());
}

TriangularKernel_RWG_SS::~TriangularKernel_RWG_SS() {
}

//virtual
dcomplex TriangularKernel_RWG_SS::rwg_value_() const noexcept {

    // The cross_g can be optimized but this does not make sense.
    double cross_g[3];
    vector_cross(Rpq_, cp_g_j_, cross_g);
    return vector_dot(cp_f_i_, cross_g) * lp_i_ * lq_j_ / 12.0;
}

///////////////////////////////////////////////////////////////////////////////

TriangularKernel_nxRWG_SS::TriangularKernel_nxRWG_SS():
TriangularKernel_RWG(),
np_{0.0, 0.0, 0.0} {
    up_green_func_.reset(new GradHelmgolzGreenFunc());
}

TriangularKernel_nxRWG_SS::~TriangularKernel_nxRWG_SS() {
}

//virtual
void TriangularKernel_nxRWG_SS::precompute_vertex_dependent_data_() noexcept {

    TriangularKernel_RWG::precompute_vertex_dependent_data_();

    const double  ev_1p[3] = {rp3_[0] - rp2_[0], rp3_[1] - rp2_[1], rp3_[2] - rp2_[2]};
    const double  ev_2p[3] = {rp1_[0] - rp3_[0], rp1_[1] - rp3_[1], rp1_[2] - rp3_[2]};

    vector_cross(ev_1p, ev_2p, np_);
    const double t_norm_np = sqrt(vector_dot(np_, np_));
    np_[0] /= t_norm_np;
    np_[1] /= t_norm_np;
    np_[2] /= t_norm_np;
}

//virtual
dcomplex TriangularKernel_nxRWG_SS::rwg_value_() const noexcept {

    double cross_f[3], cross_g[3];
    vector_cross(np_, cp_f_i_, cross_f);
    vector_cross(Rpq_, cp_g_j_, cross_g);
    return vector_dot(cross_f, cross_g) * lp_i_ * lq_j_ / 12.0;
}

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file

