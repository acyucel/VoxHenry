#ifndef _DIRECTFN_QUADRILATERAL_VOXHENRY_H_
#define _DIRECTFN_QUADRILATERAL_VOXHENRY_H_

#include <memory>
#include "directfn_defs.h"
#include "directfn_common.h"
#include "directfn_interface.h"
#include "directfn_contour.h"
#include "directfn_algorithm_voxhenry.h"
#include "directfn_kernel_quad_voxhenry.h"

using  std::unique_ptr;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_st_plan_voxhenry(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, const double nq[3], const double np[3],
                          const int lp, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_ea_plan_voxhenry(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, const double nq[3], const double np[3],
                          const int lp, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_va_plan_voxhenry(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const double r7[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, const double nq[3], const double np[3],
                          const int lp, dcomplex * const cp_data);


///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  // _DIRECTFN_QUADRILATERAL_VOXHENRY_H_

// End of the file
