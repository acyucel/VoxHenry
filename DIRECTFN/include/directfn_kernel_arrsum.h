#ifndef _DIRECTFN_KERNEL_ARRAY_H_
#define _DIRECTFN_KERNEL_ARRAY_H_

#include <memory>
#include <vector>
#include "directfn_common.h"

using  std::vector;
using  std::shared_ptr;
using  std::unique_ptr;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class KernelArrayInterface {
public:
    KernelArrayInterface();
    virtual ~KernelArrayInterface();

    KernelArrayInterface(const KernelArrayInterface & ) = delete;
    KernelArrayInterface(KernelArrayInterface && ) = delete;
    KernelArrayInterface & operator = (const KernelArrayInterface & ) = delete;
    KernelArrayInterface & operator = (KernelArrayInterface && ) = delete;

    void  setup(ParticularKernel * const ) noexcept;

    virtual void  assign_Ipsi_1_to(dcomplex * const Isd_k) const noexcept = 0;

    /*! Set internal accumulators (1,4,9,16) to zero value */
    virtual void  nullify_Ipsi_1() noexcept = 0;
    virtual void  nullify_Ieta_2() noexcept = 0;
    virtual void  nullify_Ilam_3() noexcept = 0;
    virtual void  nullify_Irho_4() noexcept = 0;

    /*! Add (accumulate) values into internal buffer */
    virtual void  accumulate_Ipsi_1(const double dfactor) noexcept = 0;
    virtual void  accumulate_Ieta_2(const double dfactor) noexcept = 0;
    virtual void  accumulate_Ilam_3(const double dfactor) noexcept = 0;
    virtual void  accumulate_Irho_4(const double dfactor) noexcept = 0;

    virtual void  multiply_Ipsi_1(const double mlt_dfactor) noexcept = 0;
    virtual void  multiply_Ieta_2(const double mlt_dfactor) noexcept = 0;
    virtual void  multiply_Ilam_3(const double mlt_dfactor) noexcept = 0;
    virtual void  multiply_Irho_4(const double mlt_dfactor) noexcept = 0;

protected:
    ParticularKernel * pkernel_;
};

///////////////////////////////////////////////////////////////////////////////

/*! In case Kernel VA, EA contains zero elements, apply inheritance */
template <typename ParticularKernel>
class KernelArray : public KernelArrayInterface<ParticularKernel> {
public:
    KernelArray();
    virtual ~KernelArray();

    virtual void  assign_Ipsi_1_to(dcomplex * const Isd_k) const noexcept;

    virtual void  nullify_Ipsi_1() noexcept;
    virtual void  nullify_Ieta_2() noexcept;
    virtual void  nullify_Ilam_3() noexcept;
    virtual void  nullify_Irho_4() noexcept;

    /*! Add (accumulate) values into internal buffer */
    virtual void  accumulate_Ipsi_1(const double add_dfactor) noexcept;
    virtual void  accumulate_Ieta_2(const double add_dfactor) noexcept;
    virtual void  accumulate_Ilam_3(const double add_dfactor) noexcept;
    virtual void  accumulate_Irho_4(const double add_dfactor) noexcept;

    virtual void  multiply_Ipsi_1(const double mlt_dfactor) noexcept;
    virtual void  multiply_Ieta_2(const double mlt_dfactor) noexcept;
    virtual void  multiply_Ilam_3(const double mlt_dfactor) noexcept;
    virtual void  multiply_Irho_4(const double mlt_dfactor) noexcept;

protected:
    /*! Vectorized arrays for integration of i-th kernels. */
    dcomplex I_psi_1_[ParticularKernel::constexpr_size()];
    dcomplex I_eta_2_[ParticularKernel::constexpr_size()];
    dcomplex I_lam_3_[ParticularKernel::constexpr_size()];
    dcomplex I_rho_4_[ParticularKernel::constexpr_size()];

    /*! Index jumper for those kernels which are exactly zero.
     *  0 for VA-SS, 1,3 for EA-SS. Thay are not summed at all*/
    size_t jmp_idx_[ParticularKernel::constexpr_size()];

    /*! Size of the jmp_idx array */
    size_t reduced_idx_sz_;

    void  index_reduction_(const std::vector<size_t> jmp_vec) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class KernelArray_TR_VA_SS final : public KernelArray<ParticularKernel> {
public:
    KernelArray_TR_VA_SS();
    virtual ~KernelArray_TR_VA_SS();
};

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class KernelArray_TR_EA_SS final : public KernelArray<ParticularKernel> {
public:
    KernelArray_TR_EA_SS();
    virtual ~KernelArray_TR_EA_SS();
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif    // _DIRECTFN_KERNEL_ARRAY_H_

// End of the file





