#ifndef _DIRECTFN_KERNEL_QUAD_GEOM_H_
#define _DIRECTFN_KERNEL_QUAD_GEOM_H_

#include <memory>
#include "directfn_kernel_base.h"

using  std::shared_ptr;
using  std::unique_ptr;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralPlanarKernel : public virtual QuadrilateralKernel {
public:
    QuadrilateralPlanarKernel();
    virtual ~QuadrilateralPlanarKernel();

    QuadrilateralPlanarKernel(const QuadrilateralPlanarKernel & ) = delete;
    QuadrilateralPlanarKernel(QuadrilateralPlanarKernel && ) = delete;
    QuadrilateralPlanarKernel & operator = (const QuadrilateralPlanarKernel & ) = delete;
    QuadrilateralPlanarKernel & operator = (QuadrilateralPlanarKernel & ) = delete;

    virtual void update_rp(const double uvxi_p[3]) noexcept;
    virtual void update_rq(const double uvxi_q[3]) noexcept;

    /*! Setup vertexes for ST */
    virtual void  set(const SingularContour3xn &  contour_xpts) noexcept;

    virtual void debug_print() const noexcept;

protected:
    /*! Vertices coordinates for the first quadrilateral */
    double rp1_[3], rp2_[3], rp3_[3], rp4_[3];

    /*! Vertices coordinates for the second quadrilateral */
    double rq1_[3], rq2_[3], rq3_[3], rq4_[3];

    /*! Linear compbinations of rp1_.. rp4_ with constant coefficients . */
    double  rp_pppp_[3];
    double  rp_mppm_[3];
    double  rp_mmpp_[3];
    double  rp_pmpm_[3];

    /*! Linear compbinations of rq1_.. rq4_ with constant coefficients . */
    double  rq_pppp_[3];
    double  rq_mppm_[3];
    double  rq_mmpp_[3];
    double  rq_pmpm_[3];

    /*! rp/rq: rp1_ + rp2_ + rp3_ + rp4_, etc...  */
    virtual void   precompute_vertex_dependent_data_() noexcept;

    /*! Derivatives depend on vertex coordinates and up, vp parameters (of transformation) */
    virtual void   precompute_rp_derivatives_() noexcept;
    /*! Derivatives depend on vertex coordinates and uq, vq parameters (of transformation) */
    virtual void   precompute_rq_derivatives_() noexcept;

    void  set_4_pts_(const SingularContour3xn &  contour_xpts) noexcept;
    void  set_6_pts_(const SingularContour3xn &  contour_xpts) noexcept;
    void  set_7_pts_(const SingularContour3xn &  contour_xpts) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralCurvilinearKernel : public virtual QuadrilateralKernel {
public:
    QuadrilateralCurvilinearKernel();
    virtual ~QuadrilateralCurvilinearKernel();

    QuadrilateralCurvilinearKernel(const QuadrilateralCurvilinearKernel & ) = delete;
    QuadrilateralCurvilinearKernel(QuadrilateralCurvilinearKernel && ) = delete;
    QuadrilateralCurvilinearKernel & operator = (const QuadrilateralCurvilinearKernel & ) = delete;
    QuadrilateralCurvilinearKernel & operator = (QuadrilateralCurvilinearKernel & ) = delete;

    virtual void update_rp(const double uvxi_p[3]) noexcept;
    virtual void update_rq(const double uvxi_q[3]) noexcept;

    /*! Setup vertexes for ST */
    virtual void  set(const SingularContour3xn &  contour_xpts) noexcept;

    virtual void debug_print() const noexcept;

protected:

    /*! Vertices coordinates for the first quadrilateral */
    double rp11_[3], rp21_[3], rp31_[3];
    double rp12_[3], rp22_[3], rp32_[3];
    double rp13_[3], rp23_[3], rp33_[3];

    /*! Vertices coordinates for the second quadrilateral */
    double rq11_[3], rq21_[3], rq31_[3];
    double rq12_[3], rq22_[3], rq32_[3];
    double rq13_[3], rq23_[3], rq33_[3];

    /*! rp/rq: rp1_ + rp2_ + rp3_ + rp4_, etc...  */
    virtual void   precompute_vertex_dependent_data_() noexcept;

    /*! Derivatives depend on vertex coordinates and up, vp parameters (of transformation) */
    virtual void   precompute_rp_derivatives_() noexcept;
    /*! Derivatives depend on vertex coordinates and uq, vq parameters (of transformation) */
    virtual void   precompute_rq_derivatives_() noexcept;

    /*! Self term  */
    void  set_9_pts_ (const SingularContour3xn &  contour_xpts) noexcept;
    void  set_15_pts_(const SingularContour3xn &  contour_xpts) noexcept;
    void  set_17_pts_(const SingularContour3xn &  contour_xpts) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}   // End of the namespace

#endif    // _DIRECTFN_KERNEL_QUAD_GEOM_H_

// End of the file




