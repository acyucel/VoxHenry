#ifndef _DIRECTFN_KERNEL_TRI_H_
#define _DIRECTFN_KERNEL_TRI_H_

#include <memory>
#include "directfn_kernel_base.h"

using  std::shared_ptr;
using  std::unique_ptr;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

class TriangularKernel : public AbstractKernel {
public:
    TriangularKernel();
    virtual ~TriangularKernel();

    TriangularKernel(const TriangularKernel & ) = delete;
    TriangularKernel(TriangularKernel && ) = delete;
    TriangularKernel & operator = (const TriangularKernel & ) = delete;
    TriangularKernel & operator = (TriangularKernel & ) = delete;

    /*! Recomputes rp, by the input uvxi_p */
    virtual void update_rp(const double uvxi_p[3]) noexcept;
    /*! Recomputes rq, by the input uvxi_q */
    virtual void update_rq(const double uvxi_q[3]) noexcept;

    /*! Setup vertexes for ST */
    virtual void  set(const SingularContour3xn &  contour_xpts) noexcept;

    virtual void debug_print() const noexcept;

protected:

    /*! Vertices coordinates for the first triangle */
    double rp1_[3], rp2_[3], rp3_[3];

    /*! Vertices coordinates for the second triangle */
    double rq1_[3], rq2_[3], rq3_[3];

    void  set_3_pts_(const SingularContour3xn &  contour_xpts) noexcept;
    void  set_4_pts_(const SingularContour3xn &  contour_xpts) noexcept;
    void  set_5_pts_(const SingularContour3xn &  contour_xpts) noexcept;

    virtual void   precompute_rp_derivatives_() noexcept;
    virtual void   precompute_rq_derivatives_() noexcept;

    virtual double calc_local_Jacobian_() noexcept;
};

/////////////////////////////////////////////////////////////////////////////////

class TriangularKernel_Constant : public TriangularKernel {
public:

    TriangularKernel_Constant();
    virtual ~TriangularKernel_Constant();

    TriangularKernel_Constant(const TriangularKernel_Constant & ) = delete;
    TriangularKernel_Constant(TriangularKernel_Constant && ) = delete;
    TriangularKernel_Constant & operator = (const TriangularKernel_Constant & ) = delete;
    TriangularKernel_Constant & operator = (TriangularKernel_Constant && ) = delete;

    virtual size_t size() const noexcept;

    /*! constexpr requires function definition in h file whereever it is used */
    static constexpr size_t constexpr_size() noexcept {return 1;}

    virtual double precomputed_jacobian() const noexcept;

protected:
    /*! Factor which is prehashed before  */
    double jacobian_;

    /*! Precomputes the Jacobian */
    virtual void precompute_vertex_dependent_data_() noexcept = 0;
    /*! Empty for Constant */
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    /*! Returns Helmgoltz Green function, index == 1 only. */
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

///----------------------------------------------------------------------------

class TriangularKernel_Constant_ST final : public TriangularKernel_Constant {
public:
    TriangularKernel_Constant_ST() = default;
    virtual ~TriangularKernel_Constant_ST() = default;
protected:
    virtual void precompute_vertex_dependent_data_() noexcept;
};

///----------------------------------------------------------------------------

class TriangularKernel_Constant_EA final : public TriangularKernel_Constant {
public:
    TriangularKernel_Constant_EA() = default;
    virtual ~TriangularKernel_Constant_EA() = default;
protected:
    virtual void precompute_vertex_dependent_data_() noexcept;
};

///----------------------------------------------------------------------------

class TriangularKernel_Constant_VA final : public TriangularKernel_Constant {
public:
    TriangularKernel_Constant_VA() = default;
    virtual ~TriangularKernel_Constant_VA() = default;
protected:
    virtual void precompute_vertex_dependent_data_() noexcept;
};

///////////////////////////////////////////////////////////////////////////////

// Forward declaration
class TriangularKernel_RWG;

using ptr_MF_RWG_VALUE = dcomplex (TriangularKernel_RWG::*)() const; // noexcept;


class TriangularKernel_RWG : public TriangularKernel {
public:

    TriangularKernel_RWG();
    virtual ~TriangularKernel_RWG();

    TriangularKernel_RWG(const TriangularKernel_RWG & ) = delete;
    TriangularKernel_RWG(TriangularKernel_RWG && ) = delete;
    TriangularKernel_RWG & operator = (const TriangularKernel_RWG & ) = delete;
    TriangularKernel_RWG & operator = (TriangularKernel_RWG && ) = delete;

    virtual size_t size() const noexcept;
    /*! The constexpr is inline. */
    static constexpr size_t constexpr_size() noexcept {return 9;}

    /*! Returns 1 */
    virtual double precomputed_jacobian() const noexcept;

protected:
    double lp_[3];
    double lq_[3];

    double f_1_[3], f_2_[3], f_3_[3];
    double g_1_[3], g_2_[3], g_3_[3];

    mutable const double * cp_f_i_;
    mutable const double * cp_g_j_;
    mutable double lp_i_;
    mutable double lq_j_;

protected:
    /*! After vertexes setup they does not change. */
    virtual void precompute_vertex_dependent_data_() noexcept;
    /*! Setup f_1_, f_2_, f_3_, g_1_, g_2_, g_3_ */
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:

    /*! Array of 1-9 function pointers. Binds dynamically
     *  input kernel index and static template parameter. */
    ptr_MF_RWG_VALUE   pmf_rwg_val_[9];

    /*! Actual kernel value */
    virtual dcomplex specific_value_(const size_t index) const noexcept;

    template <size_t ker_index>  dcomplex value_ind_() const noexcept;

    /*! WS, SS, nx_SS */
    virtual dcomplex rwg_value_() const noexcept = 0;
};

///////////////////////////////////////////////////////////////////////////////

class TriangularKernel_RWG_WS final : public TriangularKernel_RWG {
public:

    TriangularKernel_RWG_WS();
    ~TriangularKernel_RWG_WS();

    TriangularKernel_RWG_WS(const TriangularKernel_RWG_WS & ) = delete;
    TriangularKernel_RWG_WS(TriangularKernel_RWG_WS && ) = delete;
    TriangularKernel_RWG_WS & operator = (const TriangularKernel_RWG_WS & ) = delete;
    TriangularKernel_RWG_WS & operator = (TriangularKernel_RWG_WS && ) = delete;

private:
    virtual dcomplex rwg_value_() const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class TriangularKernel_RWG_SS final : public TriangularKernel_RWG {
public:
    TriangularKernel_RWG_SS();
    ~TriangularKernel_RWG_SS();

    TriangularKernel_RWG_SS(const TriangularKernel_RWG_SS & ) = delete;
    TriangularKernel_RWG_SS(TriangularKernel_RWG_SS && ) = delete;
    TriangularKernel_RWG_SS & operator = (const TriangularKernel_RWG_SS & ) = delete;
    TriangularKernel_RWG_SS & operator = (TriangularKernel_RWG_SS && ) = delete;
private:
    virtual dcomplex rwg_value_() const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class TriangularKernel_nxRWG_SS final : public TriangularKernel_RWG {
public:
    TriangularKernel_nxRWG_SS();
    ~TriangularKernel_nxRWG_SS();

    TriangularKernel_nxRWG_SS(const TriangularKernel_nxRWG_SS & ) = delete;
    TriangularKernel_nxRWG_SS(TriangularKernel_nxRWG_SS && ) = delete;
    TriangularKernel_nxRWG_SS & operator = (const TriangularKernel_nxRWG_SS & ) = delete;
    TriangularKernel_nxRWG_SS & operator = (TriangularKernel_nxRWG_SS && ) = delete;
protected:
    virtual void precompute_vertex_dependent_data_() noexcept;
private:
    double np_[3];

    virtual dcomplex rwg_value_() const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif    // _DIRECTFN_KERNEL_H_

// End of the file



