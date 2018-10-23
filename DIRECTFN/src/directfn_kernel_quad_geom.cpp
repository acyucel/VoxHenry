#include <iostream>
#include <iomanip>
#include "directfn_kernel_quad_geom.h"
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

QuadrilateralPlanarKernel::QuadrilateralPlanarKernel():
QuadrilateralKernel(),
rp1_{0.0, 0.0, 0.0},
rp2_{0.0, 0.0, 0.0},
rp3_{0.0, 0.0, 0.0},
rp4_{0.0, 0.0, 0.0},
rq1_{0.0, 0.0, 0.0},
rq2_{0.0, 0.0, 0.0},
rq3_{0.0, 0.0, 0.0},
rq4_{0.0, 0.0, 0.0},
rp_pppp_{0.0, 0.0, 0.0},
rp_mppm_{0.0, 0.0, 0.0},
rp_mmpp_{0.0, 0.0, 0.0},
rp_pmpm_{0.0, 0.0, 0.0},
rq_pppp_{0.0, 0.0, 0.0},
rq_mppm_{0.0, 0.0, 0.0},
rq_mmpp_{0.0, 0.0, 0.0},
rq_pmpm_{0.0, 0.0, 0.0} {

}

// virtual
QuadrilateralPlanarKernel::~QuadrilateralPlanarKernel() {
}

void QuadrilateralPlanarKernel::update_rp(const double uvxi_p[3]) noexcept {

    // up_  and  vp_ are hashed here to be used in the Kernel computation
    up_ = uvxi_p[0];
    vp_ = uvxi_p[1];
    const double uv_pp = up_ * vp_;

    for (size_t i = 0; i < 3; ++i) {
        rp_crnt_[i] = rp_pppp_[i] + up_ * rp_mppm_[i] + vp_ * rp_mmpp_[i] + uv_pp * rp_pmpm_[i];
    }

    precompute_rp_derivatives_();
}

void QuadrilateralPlanarKernel::update_rq(const double uvxi_q[3]) noexcept {

    // uq_  and  vq_ are hashed here to be used in the Kernel computation
    uq_ = uvxi_q[0];
    vq_ = uvxi_q[1];
    const double uv_qq = uq_ * vq_;

    for (size_t i = 0; i < 3; ++i) {
        rq_crnt_[i] = rq_pppp_[i] + uq_ * rq_mppm_[i] + vq_ * rq_mmpp_[i] + uv_qq * rq_pmpm_[i];
    }

    precompute_rq_derivatives_();
}

//virtual
void QuadrilateralPlanarKernel::set(const SingularContour3xn &  contour_xpts) noexcept {

    if (4 == contour_xpts.length()) {
        set_4_pts_(contour_xpts);
        return;
    }
    if (6 == contour_xpts.length()) {
        set_6_pts_(contour_xpts);
        return;
    }
    if (7 == contour_xpts.length()) {
        set_7_pts_(contour_xpts);
        return;
    }
    cerr << "Too few vertexes to setup quadrilateral element" << endl;
}

//virtual
void QuadrilateralPlanarKernel::precompute_rp_derivatives_() noexcept {

    for (size_t i = 0; i < 3; ++i) {
        ru_p_[i] = rp_mppm_[i] + vp_ * rp_pmpm_[i];
        rv_p_[i] = rp_mmpp_[i] + up_ * rp_pmpm_[i];
    }
}

//virtual
void QuadrilateralPlanarKernel::precompute_rq_derivatives_() noexcept {

    for (size_t i = 0; i < 3; ++i) {
        ru_q_[i] = rq_mppm_[i] + vq_ * rq_pmpm_[i];
        rv_q_[i] = rq_mmpp_[i] + uq_ * rq_pmpm_[i];
    }
}

void QuadrilateralPlanarKernel::set_4_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r1 = contour_xpts(0);
    const double * r2 = contour_xpts(1);
    const double * r3 = contour_xpts(2);
    const double * r4 = contour_xpts(3);

    for (size_t i = 0; i < 3; ++i) {

        rp1_[i] = r1[i];
        rp2_[i] = r2[i];
        rp3_[i] = r3[i];
        rp4_[i] = r4[i];

        rq1_[i] = r1[i];
        rq2_[i] = r2[i];
        rq3_[i] = r3[i];
        rq4_[i] = r4[i];
    }
    precompute_vertex_dependent_data_();
}

void QuadrilateralPlanarKernel::set_6_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r1 = contour_xpts(0);
    const double * r2 = contour_xpts(1);
    const double * r3 = contour_xpts(2);
    const double * r4 = contour_xpts(3);
    const double * r5 = contour_xpts(4);
    const double * r6 = contour_xpts(5);

    for (size_t i = 0; i < 3; ++i) {

        rp1_[i] = r4[i];
        rp2_[i] = r3[i];
        rp3_[i] = r6[i];
        rp4_[i] = r5[i];

        rq1_[i] = r3[i];
        rq2_[i] = r4[i];
        rq3_[i] = r1[i];
        rq4_[i] = r2[i];
    }
    precompute_vertex_dependent_data_();
}

void QuadrilateralPlanarKernel::set_7_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r1 = contour_xpts(0);
    const double * r2 = contour_xpts(1);
    const double * r3 = contour_xpts(2);
    const double * r4 = contour_xpts(3);
    const double * r5 = contour_xpts(4);
    const double * r6 = contour_xpts(5);
    const double * r7 = contour_xpts(6);

    for (size_t i = 0; i < 3; ++i) {

        rp1_[i] = r3[i];
        rp2_[i] = r5[i];
        rp3_[i] = r6[i];
        rp4_[i] = r7[i];

        rq1_[i] = r3[i];
        rq2_[i] = r4[i];
        rq3_[i] = r1[i];
        rq4_[i] = r2[i];
    }
    precompute_vertex_dependent_data_();
}

//virtual
void QuadrilateralPlanarKernel::debug_print() const noexcept {

    cout << "rp1: " << setw(8) << rp1_[0] << setw(8) << rp1_[1] << setw(8) << rp1_[2] << endl;
    cout << "rp2: " << setw(8) << rp2_[0] << setw(8) << rp2_[1] << setw(8) << rp2_[2] << endl;
    cout << "rp3: " << setw(8) << rp3_[0] << setw(8) << rp3_[1] << setw(8) << rp3_[2] << endl;
    cout << "rp4: " << setw(8) << rp4_[0] << setw(8) << rp4_[1] << setw(8) << rp4_[2] << endl;
    cout << endl;
    cout << "rq1: " << setw(8) << rq1_[0] << setw(8) << rq1_[1] << setw(8) << rq1_[2] << endl;
    cout << "rq2: " << setw(8) << rq2_[0] << setw(8) << rq2_[1] << setw(8) << rq2_[2] << endl;
    cout << "rq3: " << setw(8) << rq3_[0] << setw(8) << rq3_[1] << setw(8) << rq3_[2] << endl;
    cout << "rq4: " << setw(8) << rq4_[0] << setw(8) << rq4_[1] << setw(8) << rq4_[2] << endl;

    if ( nullptr != up_green_func_) {
        up_green_func_->debug_print();
    } else {
        cout << "Green function is empty..." << endl;
    }
}

//virtual
void QuadrilateralPlanarKernel::precompute_vertex_dependent_data_() noexcept {

    const double one_forth = 0.25;

    for (size_t i = 0; i < 3; ++i) {

        rp_pppp_[i] = one_forth * ( rp1_[i] + rp2_[i] + rp3_[i] + rp4_[i]);
        rp_mppm_[i] = one_forth * (-rp1_[i] + rp2_[i] + rp3_[i] - rp4_[i]);
        rp_mmpp_[i] = one_forth * (-rp1_[i] - rp2_[i] + rp3_[i] + rp4_[i]);
        rp_pmpm_[i] = one_forth * ( rp1_[i] - rp2_[i] + rp3_[i] - rp4_[i]);

        rq_pppp_[i] = one_forth * ( rq1_[i] + rq2_[i] + rq3_[i] + rq4_[i]);
        rq_mppm_[i] = one_forth * (-rq1_[i] + rq2_[i] + rq3_[i] - rq4_[i]);
        rq_mmpp_[i] = one_forth * (-rq1_[i] - rq2_[i] + rq3_[i] + rq4_[i]);
        rq_pmpm_[i] = one_forth * ( rq1_[i] - rq2_[i] + rq3_[i] - rq4_[i]);
    }
}

///////////////////////////////////////////////////////////////////////////////

QuadrilateralCurvilinearKernel::QuadrilateralCurvilinearKernel():
QuadrilateralKernel(),
rp11_{0.0, 0.0, 0.0},
rp21_{0.0, 0.0, 0.0},
rp31_{0.0, 0.0, 0.0},
rp12_{0.0, 0.0, 0.0},
rp22_{0.0, 0.0, 0.0},
rp32_{0.0, 0.0, 0.0},
rp13_{0.0, 0.0, 0.0},
rp23_{0.0, 0.0, 0.0},
rp33_{0.0, 0.0, 0.0},
rq11_{0.0, 0.0, 0.0},
rq21_{0.0, 0.0, 0.0},
rq31_{0.0, 0.0, 0.0},
rq12_{0.0, 0.0, 0.0},
rq22_{0.0, 0.0, 0.0},
rq32_{0.0, 0.0, 0.0},
rq13_{0.0, 0.0, 0.0},
rq23_{0.0, 0.0, 0.0},
rq33_{0.0, 0.0, 0.0} {

}

// virtual
QuadrilateralCurvilinearKernel::~QuadrilateralCurvilinearKernel() {

}

void QuadrilateralCurvilinearKernel::update_rp(const double uvxi_p[3]) noexcept {

    // up_  and  vp_ are hashed here to be used in the Kernel computation
    up_ = uvxi_p[0];
    vp_ = uvxi_p[1];

    const double N11 =  0.25 * (1.0 - up_) * (1.0 - vp_) * up_ * vp_;
    const double N21 = -0.5  * (1.0 - up_ * up_) * (1.0 - vp_) * vp_;
    const double N31 = -0.25 * (1.0 + up_) * (1.0 - vp_) * up_ * vp_;

    const double N12 = -0.5 * (1.0 - up_) * (1.0 - vp_ * vp_) * up_;
    const double N22 =        (1.0 - up_ * up_) * (1.0 - vp_ * vp_);
    const double N32 =  0.5 * (1.0 + up_) * (1.0 - vp_ * vp_) * up_;

    const double N13 = -0.25 * (1.0 - up_) * (1.0 + vp_) * up_ * vp_;
    const double N23 =  0.5  * (1.0 - up_ * up_) * (1.0 + vp_) * vp_;
    const double N33 =  0.25 * (1.0 + up_) * (1.0 + vp_) * up_ * vp_;

    for (size_t i = 0; i < 3; ++i) {
        rp_crnt_[i] = rp11_[i] * N11 + rp21_[i] * N21 + rp31_[i] * N31 +
                      rp12_[i] * N12 + rp22_[i] * N22 + rp32_[i] * N32 +
                      rp13_[i] * N13 + rp23_[i] * N23 + rp33_[i] * N33;
    }
    precompute_rp_derivatives_();
}

void QuadrilateralCurvilinearKernel::update_rq(const double uvxi_q[3]) noexcept {

    // uq_  and  vq_ are hashed here to be used in the Kernel computation
    uq_ = uvxi_q[0];
    vq_ = uvxi_q[1];

    const double N11 =  0.25 * (1.0 - uq_) * (1.0 - vq_) * uq_ * vq_;
    const double N21 = -0.5  * (1.0 - uq_ * uq_) * (1.0 - vq_) * vq_;
    const double N31 = -0.25 * (1.0 + uq_) * (1.0 - vq_) * uq_ * vq_;

    const double N12 = -0.5 * (1.0 - uq_) * (1.0 - vq_ * vq_) * uq_;
    const double N22 =        (1.0 - uq_ * uq_) * (1.0 - vq_ * vq_);
    const double N32 =  0.5 * (1.0 + uq_) * (1.0 - vq_ * vq_) * uq_;

    const double N13 = -0.25 * (1.0 - uq_) * (1.0 + vq_) * uq_ * vq_;
    const double N23 =  0.5  * (1.0 - uq_ * uq_) * (1.0 + vq_) * vq_;
    const double N33 =  0.25 * (1.0 + uq_) * (1.0 + vq_) * uq_ * vq_;

    for (size_t i = 0; i < 3; ++i) {
        rq_crnt_[i] = rq11_[i] * N11 + rq21_[i] * N21 + rq31_[i] * N31 +
                      rq12_[i] * N12 + rq22_[i] * N22 + rq32_[i] * N32 +
                      rq13_[i] * N13 + rq23_[i] * N23 + rq33_[i] * N33;
    }
    precompute_rq_derivatives_();
}

//virtual
void QuadrilateralCurvilinearKernel::set(const SingularContour3xn &  contour_xpts) noexcept {

    if (9 == contour_xpts.length()) {
        set_9_pts_(contour_xpts);
        return;
    }
    if (15 == contour_xpts.length()) {
        set_15_pts_(contour_xpts);
        return;
    }
    if (17 == contour_xpts.length()) {
        set_17_pts_(contour_xpts);
        return;
    }
    cerr << "Too few vertexes to setup curvilinear quadrilateral element" << endl;
}

//virtual
void QuadrilateralCurvilinearKernel::debug_print() const noexcept {

    cout << "Debug print is empty for Curvilinear...." << endl;
}

//virtual
void QuadrilateralCurvilinearKernel::precompute_rp_derivatives_() noexcept {

    const double Nu11 =  0.25 * (1.0 - 2.0 * up_) * (1.0 - vp_) * vp_;
    const double Nu21 =   up_ * (1.0 - vp_) * vp_;
    const double Nu31 = -0.25 * (1.0 + 2.0 * up_) * (1.0 - vp_) * vp_;

    const double Nu12 = -0.5 * (1.0 - 2.0 * up_) * (1.0 - vp_ * vp_);
    const double Nu22 = -2.0 * up_ * (1.0 - vp_ * vp_);
    const double Nu32 =  0.5 * (1.0 + 2.0 * up_) * (1.0 - vp_ * vp_);

    const double Nu13 = -0.25 * (1.0 - 2.0 * up_) * (1.0 + vp_) * vp_;
    const double Nu23 = -up_*(1.0 + vp_) * vp_;
    const double Nu33 =  0.25 * (1.0 + 2.0 * up_) * (1.0 + vp_) * vp_;

    for (size_t i = 0; i < 3; ++i) {

        ru_p_[i] = rp11_[i] * Nu11 + rp21_[i] * Nu21 + rp31_[i] * Nu31 +
                   rp12_[i] * Nu12 + rp22_[i] * Nu22 + rp32_[i] * Nu32 +
                   rp13_[i] * Nu13 + rp23_[i] * Nu23 + rp33_[i] * Nu33;
    }

    const double Nv11 =  0.25 * (1.0 - up_) * up_ * (1.0 - 2.0 * vp_);
    const double Nv21 = -0.5  * (1.0 - up_ * up_) * (1.0 - 2.0 * vp_);
    const double Nv31 = -0.25 * (1.0 + up_) * up_ * (1.0 - 2.0 * vp_);

    const double Nv12 = (1.0 - up_) * up_ * vp_;
    const double Nv22 = -2.0*(1 - up_ * up_) * vp_;
    const double Nv32 = -(1.0 + up_) * up_ * vp_;

    const double Nv13 = -0.25 * (1.0 - up_) * up_ * (1.0 + 2.0 * vp_);
    const double Nv23 =  0.5  * (1.0 - up_ * up_) * (1.0 + 2.0 * vp_);
    const double Nv33 =  0.25 * (1.0 + up_) * up_ * (1.0 + 2.0 * vp_);


    for (size_t i = 0; i < 3; ++i) {

        rv_p_[i] = rp11_[i] * Nv11 + rp21_[i] * Nv21 + rp31_[i] * Nv31 +
                   rp12_[i] * Nv12 + rp22_[i] * Nv22 + rp32_[i] * Nv32 +
                   rp13_[i] * Nv13 + rp23_[i] * Nv23 + rp33_[i] * Nv33;
    }
}

//virtual
void QuadrilateralCurvilinearKernel::precompute_rq_derivatives_() noexcept {

    const double Nu11 =  0.25 * (1.0 - 2.0 * uq_) * (1.0 - vq_) * vq_;
    const double Nu21 =   uq_ * (1.0 - vq_) * vq_;
    const double Nu31 = -0.25 * (1.0 + 2.0 * uq_) * (1.0 - vq_) * vq_;

    const double Nu12 = -0.5 * (1.0 - 2.0 * uq_) * (1.0 - vq_ * vq_);
    const double Nu22 = -2.0 * uq_ * (1.0 - vq_ * vq_);
    const double Nu32 =  0.5 * (1.0 + 2.0 * uq_) * (1.0 - vq_ * vq_);

    const double Nu13 = -0.25 * (1.0 - 2.0 * uq_) * (1.0 + vq_) * vq_;
    const double Nu23 = -uq_*(1.0 + vq_) * vq_;
    const double Nu33 =  0.25 * (1.0 + 2.0 * uq_) * (1.0 + vq_) * vq_;

    for (size_t i = 0; i < 3; ++i) {

        ru_q_[i] = rq11_[i] * Nu11 + rq21_[i] * Nu21 + rq31_[i] * Nu31 +
                   rq12_[i] * Nu12 + rq22_[i] * Nu22 + rq32_[i] * Nu32 +
                   rq13_[i] * Nu13 + rq23_[i] * Nu23 + rq33_[i] * Nu33;
    }

    const double Nv11 =  0.25 * (1.0 - uq_) * uq_ * (1.0 - 2.0 * vq_);
    const double Nv21 = -0.5  * (1.0 - uq_ * uq_) * (1.0 - 2.0 * vq_);
    const double Nv31 = -0.25 * (1.0 + uq_) * uq_ * (1.0 - 2.0 * vq_);

    const double Nv12 = (1.0 - uq_) * uq_ * vq_;
    const double Nv22 = -2.0*(1 - uq_ * uq_) * vq_;
    const double Nv32 = -(1.0 + uq_) * uq_ * vq_;

    const double Nv13 = -0.25 * (1.0 - uq_) * uq_ * (1.0 + 2.0 * vq_);
    const double Nv23 =  0.5  * (1.0 - uq_ * uq_) * (1.0 + 2.0 * vq_);
    const double Nv33 =  0.25 * (1.0 + uq_) * uq_ * (1.0 + 2.0 * vq_);


    for (size_t i = 0; i < 3; ++i) {

        rv_q_[i] = rq11_[i] * Nv11 + rq21_[i] * Nv21 + rq31_[i] * Nv31 +
                   rq12_[i] * Nv12 + rq22_[i] * Nv22 + rq32_[i] * Nv32 +
                   rq13_[i] * Nv13 + rq23_[i] * Nv23 + rq33_[i] * Nv33;
    }
}

void QuadrilateralCurvilinearKernel::set_9_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r1 = contour_xpts(0);
    const double * r2 = contour_xpts(1);
    const double * r3 = contour_xpts(2);
    const double * r4 = contour_xpts(3);
    const double * r5 = contour_xpts(4);
    const double * r6 = contour_xpts(5);
    const double * r7 = contour_xpts(6);
    const double * r8 = contour_xpts(7);
    const double * r9 = contour_xpts(8);

    for (size_t i = 0; i < 3; ++i) {

        rp11_[i] = r1[i];
        rp21_[i] = r2[i];
        rp31_[i] = r3[i];
        rp12_[i] = r4[i];
        rp22_[i] = r5[i];
        rp32_[i] = r6[i];
        rp13_[i] = r7[i];
        rp23_[i] = r8[i];
        rp33_[i] = r9[i];

        rq11_[i] = r1[i];
        rq21_[i] = r2[i];
        rq31_[i] = r3[i];
        rq12_[i] = r4[i];
        rq22_[i] = r5[i];
        rq32_[i] = r6[i];
        rq13_[i] = r7[i];
        rq23_[i] = r8[i];
        rq33_[i] = r9[i];
    }
}

void QuadrilateralCurvilinearKernel::set_15_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r11p = contour_xpts(0);
    const double * r21p = contour_xpts(1);
    const double * r31p = contour_xpts(2);
    const double * r12p = contour_xpts(3);
    const double * r22p = contour_xpts(4);
    const double * r32p = contour_xpts(5);
    const double * r13p = contour_xpts(6);
    const double * r23p = contour_xpts(7);
    const double * r33p = contour_xpts(8);

    const double * r12q = contour_xpts(9);
    const double * r22q = contour_xpts(10);
    const double * r32q = contour_xpts(11);
    const double * r13q = contour_xpts(12);
    const double * r23q = contour_xpts(13);
    const double * r33q = contour_xpts(14);

    for (size_t i = 0; i < 3; ++i) {

        rp11_[i] = r33p[i];
        rp21_[i] = r23p[i];
        rp31_[i] = r13p[i];
        rp12_[i] = r32p[i];
        rp22_[i] = r22p[i];
        rp32_[i] = r12p[i];
        rp13_[i] = r31p[i];
        rp23_[i] = r21p[i];
        rp33_[i] = r11p[i];

        rq11_[i] = r13p[i];
        rq21_[i] = r23p[i];
        rq31_[i] = r33p[i];
        rq12_[i] = r12q[i];
        rq22_[i] = r22q[i];
        rq32_[i] = r32q[i];
        rq13_[i] = r13q[i];
        rq23_[i] = r23q[i];
        rq33_[i] = r33q[i];
    }
}

void QuadrilateralCurvilinearKernel::set_17_pts_(const SingularContour3xn &  contour_xpts) noexcept {

    const double * r11p = contour_xpts(0);
    const double * r21p = contour_xpts(1);
    const double * r31p = contour_xpts(2);
    const double * r12p = contour_xpts(3);
    const double * r22p = contour_xpts(4);
    const double * r32p = contour_xpts(5);
    const double * r13p = contour_xpts(6);
    const double * r23p = contour_xpts(7);
    const double * r33p = contour_xpts(8);

    const double * r21q = contour_xpts(9);
    const double * r31q = contour_xpts(10);
    const double * r12q = contour_xpts(11);
    const double * r22q = contour_xpts(12);
    const double * r32q = contour_xpts(13);
    const double * r13q = contour_xpts(14);
    const double * r23q = contour_xpts(15);
    const double * r33q = contour_xpts(16);

    for (size_t i = 0; i < 3; ++i) {

        rp11_[i] = r33p[i];
        rp21_[i] = r23p[i];
        rp31_[i] = r13p[i];
        rp12_[i] = r32p[i];
        rp22_[i] = r22p[i];
        rp32_[i] = r12p[i];
        rp13_[i] = r31p[i];
        rp23_[i] = r21p[i];
        rp33_[i] = r11p[i];

        rq11_[i] = r33p[i];
        rq21_[i] = r21q[i];
        rq31_[i] = r31q[i];
        rq12_[i] = r12q[i];
        rq22_[i] = r22q[i];
        rq32_[i] = r32q[i];
        rq13_[i] = r13q[i];
        rq23_[i] = r23q[i];
        rq33_[i] = r33q[i];
    }
}


//virtual
void QuadrilateralCurvilinearKernel::precompute_vertex_dependent_data_() noexcept {
}

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file
