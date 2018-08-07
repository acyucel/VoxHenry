#ifndef _DIRECTFN_KERNEL_QUAD_SCALAR_H_
#define _DIRECTFN_KERNEL_QUAD_SCALAR_H_

#include "directfn_kernel_base.h"
#include "directfn_kernel_quad_geom.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_Scalar : public virtual QuadrilateralKernel {
public:
    QuadrilateralKernel_Scalar();
    virtual ~QuadrilateralKernel_Scalar();

    QuadrilateralKernel_Scalar(const QuadrilateralKernel_Scalar & ) = delete;
    QuadrilateralKernel_Scalar(QuadrilateralKernel_Scalar && ) = delete;
    QuadrilateralKernel_Scalar & operator = (const QuadrilateralKernel_Scalar & ) = delete;
    QuadrilateralKernel_Scalar & operator = (QuadrilateralKernel_Scalar && ) = delete;

    virtual size_t size() const noexcept;
    static constexpr size_t constexpr_size() noexcept {return 1;}

    /*! Is used after N1...N4 loop computations and needed if the Jacobian is constant.
     *  Returns 1 because it varies inside the kernel. */
    virtual double precomputed_jacobian() const noexcept;

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    /*! It is computed once rp_, rq_ defined. */
    double jacobian_;

    /*! Actual kernel value */
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_PlanarScalar final : public QuadrilateralPlanarKernel,
                                               public QuadrilateralKernel_Scalar {
public:
    QuadrilateralKernel_PlanarScalar();
    virtual ~QuadrilateralKernel_PlanarScalar();

    QuadrilateralKernel_PlanarScalar(const QuadrilateralKernel_PlanarScalar & ) = delete;
    QuadrilateralKernel_PlanarScalar(QuadrilateralKernel_PlanarScalar && ) = delete;

    QuadrilateralKernel_PlanarScalar & operator = (const QuadrilateralKernel_PlanarScalar & ) = delete;
    QuadrilateralKernel_PlanarScalar & operator = (QuadrilateralKernel_PlanarScalar && ) = delete;
private:
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_CurvilinearScalar final : public QuadrilateralCurvilinearKernel,
                                                    public QuadrilateralKernel_Scalar {
public:
    QuadrilateralKernel_CurvilinearScalar();
    virtual ~QuadrilateralKernel_CurvilinearScalar();

    QuadrilateralKernel_CurvilinearScalar(const QuadrilateralKernel_CurvilinearScalar & ) = delete;
    QuadrilateralKernel_CurvilinearScalar(QuadrilateralKernel_CurvilinearScalar && ) = delete;

    QuadrilateralKernel_CurvilinearScalar & operator = (const QuadrilateralKernel_CurvilinearScalar & ) = delete;
    QuadrilateralKernel_CurvilinearScalar & operator = (QuadrilateralKernel_CurvilinearScalar && ) = delete;
private:
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  //  _DIRECTFN_KERNEL_QUAD_SCALAR_H_

// End of the file
