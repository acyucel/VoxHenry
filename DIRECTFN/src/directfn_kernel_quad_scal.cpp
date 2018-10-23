#include "directfn_kernel_quad_scal.h"
#include "directfn_greenfunc.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_Scalar::QuadrilateralKernel_Scalar():
QuadrilateralKernel(),
jacobian_(0.0) {
    up_green_func_.reset(new HelmgolzGreenFunc());
}

//virtual
QuadrilateralKernel_Scalar::~QuadrilateralKernel_Scalar() {
}

//virtual
size_t QuadrilateralKernel_Scalar::size() const noexcept {
    return QuadrilateralKernel_Scalar::constexpr_size();
}

//virtual
double QuadrilateralKernel_Scalar::precomputed_jacobian() const noexcept {
    // This is not jacobian_ used in the loop where rp_ rq_ changes.
    // It must be used if and only if the J is const for all rp and rq
    return 1.0;
}

//virtual
void QuadrilateralKernel_Scalar::precompute_rp_rq_dependent_data_()  noexcept {
    jacobian_ = calc_Jacobian_();
}


//virtual
dcomplex QuadrilateralKernel_Scalar::specific_value_(const size_t ) const noexcept {
    return  up_green_func_->value() * jacobian_;
}

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_PlanarScalar::QuadrilateralKernel_PlanarScalar():
QuadrilateralKernel(),
QuadrilateralPlanarKernel(),
QuadrilateralKernel_Scalar() {

}

//virtual
QuadrilateralKernel_PlanarScalar::~QuadrilateralKernel_PlanarScalar() {

}

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_CurvilinearScalar::QuadrilateralKernel_CurvilinearScalar():
QuadrilateralKernel(),
QuadrilateralCurvilinearKernel(),
QuadrilateralKernel_Scalar() {
}

//virtual
QuadrilateralKernel_CurvilinearScalar::~QuadrilateralKernel_CurvilinearScalar() {
}

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file


