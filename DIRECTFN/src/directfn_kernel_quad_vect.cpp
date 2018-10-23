#include "directfn_kernel_quad_vect.h"
#include "directfn_greenfunc.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<0>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t1_;
    cp_b_j_ = b1_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<1>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t1_;
    cp_b_j_ = b2_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<2>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t1_;
    cp_b_j_ = b3_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<3>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t1_;
    cp_b_j_ = b4_;
    return rooftop_value_();
}


template <> dcomplex QuadrilateralKernel_Vector::value_ind_<4>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t2_;
    cp_b_j_ = b1_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<5>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t2_;
    cp_b_j_ = b2_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<6>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t2_;
    cp_b_j_ = b3_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<7>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t2_;
    cp_b_j_ = b4_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<8>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t3_;
    cp_b_j_ = b1_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<9>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t3_;
    cp_b_j_ = b2_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<10>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t3_;
    cp_b_j_ = b3_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<11>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t3_;
    cp_b_j_ = b4_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<12>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t4_;
    cp_b_j_ = b1_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<13>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t4_;
    cp_b_j_ = b2_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<14>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t4_;
    cp_b_j_ = b3_;
    return rooftop_value_();
}

template <> dcomplex QuadrilateralKernel_Vector::value_ind_<15>() const noexcept {
    // Get corresponding addresses
    cp_t_i_ = t4_;
    cp_b_j_ = b4_;
    return rooftop_value_();
}

QuadrilateralKernel_Vector::QuadrilateralKernel_Vector():
QuadrilateralKernel(),
t1_{0.0, 0.0, 0.0},
t2_{0.0, 0.0, 0.0},
t3_{0.0, 0.0, 0.0},
t4_{0.0, 0.0, 0.0},
b1_{0.0, 0.0, 0.0},
b2_{0.0, 0.0, 0.0},
b3_{0.0, 0.0, 0.0},
b4_{0.0, 0.0, 0.0},
cp_t_i_(nullptr),
cp_b_j_(nullptr),
pmf_rt4_val_{nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
             nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr} {

    pmf_rt4_val_[0]  = &QuadrilateralKernel_Vector::value_ind_<0>;
    pmf_rt4_val_[1]  = &QuadrilateralKernel_Vector::value_ind_<1>;
    pmf_rt4_val_[2]  = &QuadrilateralKernel_Vector::value_ind_<2>;
    pmf_rt4_val_[3]  = &QuadrilateralKernel_Vector::value_ind_<3>;
    pmf_rt4_val_[4]  = &QuadrilateralKernel_Vector::value_ind_<4>;
    pmf_rt4_val_[5]  = &QuadrilateralKernel_Vector::value_ind_<5>;
    pmf_rt4_val_[6]  = &QuadrilateralKernel_Vector::value_ind_<6>;
    pmf_rt4_val_[7]  = &QuadrilateralKernel_Vector::value_ind_<7>;
    pmf_rt4_val_[8]  = &QuadrilateralKernel_Vector::value_ind_<8>;
    pmf_rt4_val_[9]  = &QuadrilateralKernel_Vector::value_ind_<9>;
    pmf_rt4_val_[10] = &QuadrilateralKernel_Vector::value_ind_<10>;
    pmf_rt4_val_[11] = &QuadrilateralKernel_Vector::value_ind_<11>;
    pmf_rt4_val_[12] = &QuadrilateralKernel_Vector::value_ind_<12>;
    pmf_rt4_val_[13] = &QuadrilateralKernel_Vector::value_ind_<13>;
    pmf_rt4_val_[14] = &QuadrilateralKernel_Vector::value_ind_<14>;
    pmf_rt4_val_[15] = &QuadrilateralKernel_Vector::value_ind_<15>;
}

//virtual
QuadrilateralKernel_Vector::~QuadrilateralKernel_Vector() {
}

////virtual
//void QuadrilateralKernel_Vector::update_rp(const double uvxi_p[3]) noexcept {

//    QuadrilateralPlanarKernel::update_rp(uvxi_p);
//    precompute_rp_only_data_();
//}

////virtual
//void QuadrilateralKernel_Vector::update_rq(const double uvxi_q[3]) noexcept {

//    QuadrilateralPlanarKernel::update_rq(uvxi_q);
//    precompute_rq_only_data_();
//}

//virtual
size_t QuadrilateralKernel_Vector::size() const noexcept {
    return constexpr_size();
}

//virtual
double QuadrilateralKernel_Vector::precomputed_jacobian() const noexcept {
    return 1.0;
}

//virtual
void QuadrilateralKernel_Vector::precompute_rp_rq_dependent_data_() noexcept {

    precompute_rp_only_data_();
    precompute_rq_only_data_();
}

void QuadrilateralKernel_Vector::precompute_rp_only_data_() noexcept {

    for (size_t i = 0; i < 3; ++i) {

        t1_[i] = (1.0 + up_) * ru_p_[i];
        t2_[i] = (1.0 - up_) * ru_p_[i];
        t3_[i] = (1.0 + vp_) * rv_p_[i];
        t4_[i] = (1.0 - vp_) * rv_p_[i];
    }
}

void QuadrilateralKernel_Vector::precompute_rq_only_data_() noexcept {

    for (size_t i = 0; i < 3; ++i) {

        b1_[i] = (1.0 + uq_) * ru_q_[i];
        b2_[i] = (1.0 - uq_) * ru_q_[i];
        b3_[i] = (1.0 + vq_) * rv_q_[i];
        b4_[i] = (1.0 - vq_) * rv_q_[i];
    }
}

//virtual
dcomplex QuadrilateralKernel_Vector::specific_value_(const size_t index) const noexcept {
    return (this->*pmf_rt4_val_[index])() * up_green_func_->value();
}

/////////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_PlanarVectorWS::QuadrilateralKernel_PlanarVectorWS():
QuadrilateralKernel(),
QuadrilateralPlanarKernel(),
QuadrilateralKernel_Vector() {
    up_green_func_.reset(new HelmgolzGreenFunc());
}

//virtual
QuadrilateralKernel_PlanarVectorWS::~QuadrilateralKernel_PlanarVectorWS() {
}

//virtual
dcomplex QuadrilateralKernel_PlanarVectorWS::rooftop_value_() const noexcept {
    return vector_dot(cp_t_i_, cp_b_j_);
}

/////////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_PlanarVectorSS::QuadrilateralKernel_PlanarVectorSS():
QuadrilateralKernel(),
QuadrilateralPlanarKernel(),
QuadrilateralKernel_Vector() {
    up_green_func_.reset(new GradHelmgolzGreenFunc());
}

//virtual
QuadrilateralKernel_PlanarVectorSS::~QuadrilateralKernel_PlanarVectorSS() {
}

//virtual
dcomplex QuadrilateralKernel_PlanarVectorSS::rooftop_value_() const noexcept {

    double cross_g[3];
    vector_cross(Rpq_, cp_b_j_, cross_g);
    return vector_dot(cp_t_i_, cross_g);
}

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_CurvilinearVectorWS::QuadrilateralKernel_CurvilinearVectorWS():
QuadrilateralKernel(),
QuadrilateralCurvilinearKernel(),
QuadrilateralKernel_Vector() {
    up_green_func_.reset(new HelmgolzGreenFunc());
}

//virtual
QuadrilateralKernel_CurvilinearVectorWS::~QuadrilateralKernel_CurvilinearVectorWS() {
}

//virtual
dcomplex QuadrilateralKernel_CurvilinearVectorWS::rooftop_value_() const noexcept {
    return vector_dot(cp_t_i_, cp_b_j_);
}

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_CurvilinearVectorSS::QuadrilateralKernel_CurvilinearVectorSS():
QuadrilateralKernel(),
QuadrilateralCurvilinearKernel(),
QuadrilateralKernel_Vector() {
    up_green_func_.reset(new GradHelmgolzGreenFunc());
}

//virtual
QuadrilateralKernel_CurvilinearVectorSS::~QuadrilateralKernel_CurvilinearVectorSS() {
}

//virtual
dcomplex QuadrilateralKernel_CurvilinearVectorSS::rooftop_value_() const noexcept {

    double cross_g[3];
    vector_cross(Rpq_, cp_b_j_, cross_g);
    return vector_dot(cp_t_i_, cross_g);
}

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file


