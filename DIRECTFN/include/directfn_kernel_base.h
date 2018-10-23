#ifndef _DIRECTFN_KERNEL_BASE_H_
#define _DIRECTFN_KERNEL_BASE_H_

#include <memory>
#include "directfn_common.h"

using  std::shared_ptr;
using  std::unique_ptr;

namespace Directfn {

class  AbstractGreenFunc;
class  SingularContour3xn;

///////////////////////////////////////////////////////////////////////////////

class AbstractKernel {
public:
    /*! Default constructor */
    AbstractKernel();
    /*! Default destructor */
    virtual ~AbstractKernel();

    AbstractKernel(const AbstractKernel & ) = delete;
    AbstractKernel(AbstractKernel && ) = delete;
    AbstractKernel & operator = (const AbstractKernel & ) = delete;
    AbstractKernel & operator = (AbstractKernel && ) = delete;

    /*! The wave number is used in the exp(i k R)  */
    void set_wavenumber(const double k0_inp) noexcept;

    /*! The maximum number of different kernel types like 1,9 or 16. */
    virtual size_t size() const noexcept = 0;

    /*! Recomputes rp by the input array of two or three parameters,
     *  rp_ is saved then. */
    virtual void update_rp(const double uvxi_p[3]) noexcept = 0;

    /*! Recomputes rq by the input array of two or three parameters
     *  rq_ is saved then. */
    virtual void update_rq(const double uvxi_q[3]) noexcept = 0;

    /*! Calculates GF value to be used in kernels loop without recalculation. */
    void  precompute_rp_rq_data() noexcept;

    /*! NVI idiom: avoid the default values for virtual functions. Item 37, [SM-55] */
    dcomplex value(const size_t index = 0) const noexcept;

    /*! The Jacobian factor which depands only on vertex parameters.
     *  In fact it is equal to 1 for all kernels excluding the Constant basis case.
     *  For RWG_xx it is equal to 1 and for Quadrilateral elements the J factor
     *  is considered inside the kernel computation. */
    virtual double precomputed_jacobian() const noexcept = 0;

    /*! The radius vector setup depends on the case of triangles, quadrangles, st, ea, va.
     *  Formal definitions used in inherited classes.*/
    virtual void  set(const SingularContour3xn &  contour_xpts) noexcept = 0;

    virtual void debug_print() const noexcept = 0;

protected:

    /*! The green function is initialized in the inherited constructors */
    unique_ptr<AbstractGreenFunc>      up_green_func_;

    /*! Parameters of integration rp and rq for Green func */
    double rp_crnt_[3], rq_crnt_[3];

    /*! Rpq_[i] = rp_[i] - rq_[i] */
    double Rpq_[3];

    /*! lp[3], lq[3], np_[3] if needed */
    virtual void precompute_vertex_dependent_data_() noexcept = 0;
    /*! f_1,2,3, g_1,2,3 */
    virtual void precompute_rp_rq_dependent_data_()  noexcept = 0;

    /*! Derivatives depend on vertex coordinates and up, vp parameters (of transformation) */
    virtual void   precompute_rp_derivatives_() noexcept = 0;
    /*! Derivatives depend on vertex coordinates and uq, vq parameters (of transformation) */
    virtual void   precompute_rq_derivatives_() noexcept = 0;

    /*! Normales are needed for nx basis functions. */
    //virtual void   calc_normales_(); // will be implemented later

private:
    /*! Reimplemented in inherited classes */
    virtual dcomplex specific_value_(const size_t index) const noexcept = 0;
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel : public AbstractKernel {
public:

    QuadrilateralKernel();
    virtual ~QuadrilateralKernel();

protected:
    /*! Parameter of transformation */
    double up_, vp_;
    double uq_, vq_;

    /*! Precomputed derivatives */
    double  ru_p_[3],  rv_p_[3];
    double  ru_q_[3],  rv_q_[3];

    double calc_Jacobian_() const noexcept;

    /*! Derivatives depend on vertex coordinates and up, vp, uq, vq parameters (of transformation) */
    /*! Derivatives depend on vertex coordinates and up, vp parameters (of transformation) */
    virtual void   precompute_rp_derivatives_() noexcept = 0;
    /*! Derivatives depend on vertex coordinates and uq, vq parameters (of transformation) */
    virtual void   precompute_rq_derivatives_() noexcept = 0;

};

///////////////////////////////////////////////////////////////////////////////

// Triangular

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif    // _DIRECTFN_KERNEL_H_

// End of the file





