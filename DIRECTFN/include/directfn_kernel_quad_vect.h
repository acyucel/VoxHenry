#ifndef _DIRECTFN_KERNEL_QUAD_VECTOR_H_
#define _DIRECTFN_KERNEL_QUAD_VECTOR_H_

#include "directfn_kernel_base.h"
#include "directfn_kernel_quad_geom.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

// Forward declaration
class QuadrilateralKernel_Vector;

using ptr_MF_RT4_VALUE = dcomplex (QuadrilateralKernel_Vector::*)() const;

class QuadrilateralKernel_Vector : public virtual QuadrilateralKernel {
public:

    QuadrilateralKernel_Vector();
    virtual ~QuadrilateralKernel_Vector();

    QuadrilateralKernel_Vector(const QuadrilateralKernel_Vector & ) = delete;
    QuadrilateralKernel_Vector(QuadrilateralKernel_Vector && ) = delete;
    QuadrilateralKernel_Vector & operator = (const QuadrilateralKernel_Vector & ) = delete;
    QuadrilateralKernel_Vector & operator = (QuadrilateralKernel_Vector && ) = delete;

    virtual size_t size() const noexcept;
    /*! Inline constexpr. */
    static constexpr size_t constexpr_size() noexcept {return 16;}

    /*! Is used after all  */
    virtual double precomputed_jacobian() const noexcept;  // 1

protected:
    double  t1_[3], t2_[3], t3_[3], t4_[4];
    double  b1_[3], b2_[3], b3_[3], b4_[4];

    mutable const double * cp_t_i_;
    mutable const double * cp_b_j_;

    /*! Setup all ru_p_,.., rv_q_, t1_,...b4_ */
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

    void precompute_rp_only_data_() noexcept;
    void precompute_rq_only_data_() noexcept;

private:
    /*! Array of 1-16 function pointers. Binds dynamically
     *  input kernel index and static template parameter. */
    ptr_MF_RT4_VALUE  pmf_rt4_val_[16];

    /*! Calls function pointer to the templated function value_ind_. */
    virtual dcomplex specific_value_(const size_t index) const noexcept;

    /*! Templated specific_value_ to be placed into array.
     *  Then call rooftop_value_ with reassigned cp_t_i_, cp_b_j_ */
    template <size_t ker_index>  dcomplex value_ind_() const noexcept;

    /*! WS, SS */
    virtual dcomplex rooftop_value_() const noexcept = 0;
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_PlanarVectorWS final : public QuadrilateralPlanarKernel,
                                                 public QuadrilateralKernel_Vector {
public:
    QuadrilateralKernel_PlanarVectorWS();
    virtual ~QuadrilateralKernel_PlanarVectorWS();

    QuadrilateralKernel_PlanarVectorWS(const QuadrilateralKernel_PlanarVectorWS & ) = delete;
    QuadrilateralKernel_PlanarVectorWS(QuadrilateralKernel_PlanarVectorWS && ) = delete;
    QuadrilateralKernel_PlanarVectorWS & operator = (const QuadrilateralKernel_PlanarVectorWS & ) = delete;
    QuadrilateralKernel_PlanarVectorWS & operator = (QuadrilateralKernel_PlanarVectorWS && ) = delete;

private:
    virtual dcomplex rooftop_value_() const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_PlanarVectorSS final : public QuadrilateralPlanarKernel,
                                                 public QuadrilateralKernel_Vector {
public:
    QuadrilateralKernel_PlanarVectorSS();
    virtual ~QuadrilateralKernel_PlanarVectorSS();

    QuadrilateralKernel_PlanarVectorSS(const QuadrilateralKernel_PlanarVectorSS & ) = delete;
    QuadrilateralKernel_PlanarVectorSS(QuadrilateralKernel_PlanarVectorSS && ) = delete;
    QuadrilateralKernel_PlanarVectorSS & operator = (const QuadrilateralKernel_PlanarVectorSS & ) = delete;
    QuadrilateralKernel_PlanarVectorSS & operator = (QuadrilateralKernel_PlanarVectorSS && ) = delete;

private:
    virtual dcomplex rooftop_value_() const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_CurvilinearVectorWS final : public QuadrilateralCurvilinearKernel,
                                                      public QuadrilateralKernel_Vector {
public:
    QuadrilateralKernel_CurvilinearVectorWS();
    virtual ~QuadrilateralKernel_CurvilinearVectorWS();

    QuadrilateralKernel_CurvilinearVectorWS(const QuadrilateralKernel_CurvilinearVectorWS & ) = delete;
    QuadrilateralKernel_CurvilinearVectorWS(QuadrilateralKernel_CurvilinearVectorWS && ) = delete;
    QuadrilateralKernel_CurvilinearVectorWS & operator = (const QuadrilateralKernel_CurvilinearVectorWS & ) = delete;
    QuadrilateralKernel_CurvilinearVectorWS & operator = (QuadrilateralKernel_CurvilinearVectorWS && ) = delete;

private:
    virtual dcomplex rooftop_value_() const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_CurvilinearVectorSS final : public QuadrilateralCurvilinearKernel,
                                                      public QuadrilateralKernel_Vector {
public:
    QuadrilateralKernel_CurvilinearVectorSS();
    virtual ~QuadrilateralKernel_CurvilinearVectorSS();

    QuadrilateralKernel_CurvilinearVectorSS(const QuadrilateralKernel_CurvilinearVectorWS & ) = delete;
    QuadrilateralKernel_CurvilinearVectorSS(QuadrilateralKernel_CurvilinearVectorSS && ) = delete;
    QuadrilateralKernel_CurvilinearVectorSS & operator = (const QuadrilateralKernel_CurvilinearVectorSS & ) = delete;
    QuadrilateralKernel_CurvilinearVectorSS & operator = (QuadrilateralKernel_CurvilinearVectorSS && ) = delete;

private:
    virtual dcomplex rooftop_value_() const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  //  _DIRECTFN_KERNEL_QUAD_VECTOR_H_

// End of the file
