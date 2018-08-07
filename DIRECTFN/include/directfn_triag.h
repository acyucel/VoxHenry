#ifndef _DIRECTFN_TRIANGULAR_ALL_H_
#define _DIRECTFN_TRIANGULAR_ALL_H_

#include <memory>
#include "directfn_defs.h"
#include "directfn_common.h"
#include "directfn_contour.h"
#include "directfn_algorithm_st.h"
#include "directfn_algorithm_ea.h"
#include "directfn_algorithm_va.h"
#include "directfn_kernel_tri.h"

using  std::unique_ptr;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

/*! cp_data must be allocated before input pass */

template <typename ParticularKernel>
int directfn_tri_st_plan(const double r1[3], const double r2[3], const double r3[3],
                         const size_t N1, const size_t N2,
                         const size_t N3, const size_t N4,
                         const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_tri_ea_plan(const double r1[3], const double r2[3],
                         const double r3[3], const double r4[3],
                         const size_t N1, const size_t N2,
                         const size_t N3, const size_t N4,
                         const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_tri_va_plan(const double r1[3], const double r2[3],
                         const double r3[3], const double r4[3],
                         const double r5[3],
                         const size_t N1, const size_t N2,
                         const size_t N3, const size_t N4,
                         const double k0, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  // _DIRECTFN_TRIANGULAR_ALL_H_

// End of the file





