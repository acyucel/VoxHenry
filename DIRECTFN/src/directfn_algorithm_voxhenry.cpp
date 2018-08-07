#include <iostream>
#include <iomanip>
#include "directfn_algorithm_voxhenry.h"
#include "directfn_kernel_quad_voxhenry.h"

using std::cout;
using std::endl;

namespace Directfn {

template <typename ParticularKernel>
Quadrilateral_EA_VH<ParticularKernel>::Quadrilateral_EA_VH() :
Quadrilateral_EA<ParticularKernel>() {

}

template <typename ParticularKernel>
Quadrilateral_EA_VH<ParticularKernel>::~Quadrilateral_EA_VH() {
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA_VH<ParticularKernel>::set_lp(const int lp) noexcept {
    this->up_kernel_->set_lp(lp);
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_EA_VH<ParticularKernel>::set_nq_np(const double nq[3], const double np[3]) noexcept {
    this->up_kernel_->set_nq_np(nq, np);
}

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Quadrilateral_ST_VH<ParticularKernel>::Quadrilateral_ST_VH() :
Quadrilateral_ST<ParticularKernel>() {

}

template <typename ParticularKernel>
Quadrilateral_ST_VH<ParticularKernel>::~Quadrilateral_ST_VH() {
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST_VH<ParticularKernel>::set_lp(const int lp) noexcept {
    this->up_kernel_->set_lp(lp);
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_ST_VH<ParticularKernel>::set_nq_np(const double nq[3], const double np[3]) noexcept {
    this->up_kernel_->set_nq_np(nq, np);
}

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
Quadrilateral_VA_VH<ParticularKernel>::Quadrilateral_VA_VH() :
Quadrilateral_VA<ParticularKernel>() {

}

template <typename ParticularKernel>
Quadrilateral_VA_VH<ParticularKernel>::~Quadrilateral_VA_VH() {
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_VA_VH<ParticularKernel>::set_lp(const int lp) noexcept {
    this->up_kernel_->set_lp(lp);
}

//virtual
template <typename ParticularKernel>
void Quadrilateral_VA_VH<ParticularKernel>::set_nq_np(const double nq[3], const double np[3]) noexcept {
    this->up_kernel_->set_nq_np(nq, np);
}

///////////////////////////////////////////////////////////////////////////////


// Instantiation of the Quadrilateral Constant VoxHenry Kernels

template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class Quadrilateral_EA_VH<QuadKer_PlanVH_VolKer3_KerTyp4>;


template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class Quadrilateral_ST_VH<QuadKer_PlanVH_VolKer3_KerTyp4>;


template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class Quadrilateral_VA_VH<QuadKer_PlanVH_VolKer3_KerTyp4>;


///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file

