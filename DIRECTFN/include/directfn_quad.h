#ifndef _DIRECTFN_QUADRILATERAL_ALL_H_
#define _DIRECTFN_QUADRILATERAL_ALL_H_

#include <memory>
#include "directfn_defs.h"
#include "directfn_common.h"
#include "directfn_interface.h"
#include "directfn_contour.h"
#include "directfn_algorithm_st.h"
#include "directfn_algorithm_ea.h"
#include "directfn_algorithm_va.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"

using  std::unique_ptr;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_st_plan(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_ea_plan(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_va_plan(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const double r7[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_st_curv(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const double r7[3], const double r8[3],
                          const double r9[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_ea_curv(const double r1[3],  const double r2[3],
                          const double r3[3],  const double r4[3],
                          const double r5[3],  const double r6[3],
                          const double r7[3],  const double r8[3],
                          const double r9[3],  const double r10[3],
                          const double r11[3], const double r12[3],
                          const double r13[3], const double r14[3],
                          const double r15[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_va_curv(const double r1[3],  const double r2[3],
                          const double r3[3],  const double r4[3],
                          const double r5[3],  const double r6[3],
                          const double r7[3],  const double r8[3],
                          const double r9[3],  const double r10[3],
                          const double r11[3], const double r12[3],
                          const double r13[3], const double r14[3],
                          const double r15[3], const double r16[3],
                          const double r17[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  // _DIRECTFN_QUADRILATERAL_ALL_H_

// End of the file
