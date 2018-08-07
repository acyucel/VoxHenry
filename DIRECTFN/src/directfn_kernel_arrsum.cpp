#include <iostream>
#include <iomanip>
#include <algorithm>
#include "directfn_interface.h"
#include "directfn_kernel_base.h"
#include "directfn_kernel_tri.h"
#include "directfn_kernel_quad_voxhenry.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"
#include "directfn_kernel_quad_voxhenry.h"

using  std::cout;
using  std::cerr;
using  std::endl;
using  std::setprecision;
using  std::setw;
using  std::to_string;
using  std::vector;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
KernelArrayInterface<ParticularKernel>::KernelArrayInterface():
pkernel_(nullptr) {
}

//virtual
template <typename ParticularKernel>
KernelArrayInterface<ParticularKernel>::~KernelArrayInterface() {
}

template <typename ParticularKernel>
void KernelArrayInterface<ParticularKernel>::setup(ParticularKernel * const inp_Kernel) noexcept {
    pkernel_ = inp_Kernel;
}


// Triangular elements
template class KernelArrayInterface<TriangularKernel_Constant_ST>;
template class KernelArrayInterface<TriangularKernel_Constant_EA>;
template class KernelArrayInterface<TriangularKernel_Constant_VA>;

template class KernelArrayInterface<TriangularKernel_RWG_WS>;
template class KernelArrayInterface<TriangularKernel_RWG_SS>;
template class KernelArrayInterface<TriangularKernel_nxRWG_SS>;

// Quadrilateral for Planar
template class KernelArrayInterface<QuadrilateralKernel_PlanarScalar>;
template class KernelArrayInterface<QuadrilateralKernel_PlanarVectorWS>;
template class KernelArrayInterface<QuadrilateralKernel_PlanarVectorSS>;

template class KernelArrayInterface<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class KernelArrayInterface<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class KernelArrayInterface<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class KernelArrayInterface<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class KernelArrayInterface<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class KernelArrayInterface<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class KernelArrayInterface<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class KernelArrayInterface<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class KernelArrayInterface<QuadrilateralKernel_CurvilinearScalar>;
template class KernelArrayInterface<QuadrilateralKernel_CurvilinearVectorWS>;
template class KernelArrayInterface<QuadrilateralKernel_CurvilinearVectorSS>;

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
KernelArray<ParticularKernel>::KernelArray():
KernelArrayInterface<ParticularKernel>(),
reduced_idx_sz_(ParticularKernel::constexpr_size()) {

    vector<size_t> exclude_idxs(0);
    index_reduction_(exclude_idxs);
}

template <typename ParticularKernel>
KernelArray<ParticularKernel>::~KernelArray() {
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::assign_Ipsi_1_to(dcomplex * const Isd_k) const noexcept {

    // Set unused values to zero
    for (size_t i = 0; i < ParticularKernel::constexpr_size(); ++i) {
        Isd_k[i] = dcomplex(0.0);
    }
    // Load non-zero values
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        Isd_k[jmp_idx_[i]] = I_psi_1_[i];
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::nullify_Ipsi_1() noexcept {

    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_psi_1_[i] = dcomplex(0.0, 0.0);
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::nullify_Ieta_2() noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_eta_2_[i] = dcomplex(0.0, 0.0);
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::nullify_Ilam_3() noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_lam_3_[i] = dcomplex(0.0, 0.0);
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::nullify_Irho_4() noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_rho_4_[i] = dcomplex(0.0, 0.0);
    }
}

// virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::accumulate_Ipsi_1(const double add_dfactor) noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_psi_1_[i] += add_dfactor * I_eta_2_[i];
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::accumulate_Ieta_2(const double add_dfactor) noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_eta_2_[i] += add_dfactor * I_lam_3_[i];
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::accumulate_Ilam_3(const double add_dfactor) noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_lam_3_[i] += add_dfactor * I_rho_4_[i];
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::accumulate_Irho_4(const double add_dfactor) noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        // Shifted adressation
        I_rho_4_[i] += add_dfactor * this->pkernel_->value(jmp_idx_[i]);
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::multiply_Ipsi_1(const double mlt_dfactor) noexcept {

    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_psi_1_[i] *= mlt_dfactor;
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::multiply_Ieta_2(const double mlt_dfactor) noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_eta_2_[i] *= mlt_dfactor;
    }
}

//virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::multiply_Ilam_3(const double mlt_dfactor) noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_lam_3_[i] *= mlt_dfactor;
    }
}

// virtual
template <typename ParticularKernel>
void KernelArray<ParticularKernel>::multiply_Irho_4(const double mlt_dfactor) noexcept {
    for (size_t i = 0; i < reduced_idx_sz_; ++i) {
        I_rho_4_[i] *= mlt_dfactor;
    }
}

template <typename ParticularKernel>
void KernelArray<ParticularKernel>::index_reduction_(const std::vector<size_t> jmp_vec) noexcept {

    reduced_idx_sz_ = ParticularKernel::constexpr_size() - jmp_vec.size();

    for (size_t i_std = 0, k_jmp = 0; i_std < ParticularKernel::constexpr_size(); ++i_std) {
        std::vector<size_t>::const_iterator  itrf = std::find(jmp_vec.begin(), jmp_vec.end(), i_std);
        if (itrf != jmp_vec.cend()) {
            continue;
        }
        jmp_idx_[k_jmp++] = i_std;
    }
}


// Triangular elements
template class KernelArray<TriangularKernel_Constant_ST>;
template class KernelArray<TriangularKernel_Constant_EA>;
template class KernelArray<TriangularKernel_Constant_VA>;

template class KernelArray<TriangularKernel_RWG_WS>;
template class KernelArray<TriangularKernel_RWG_SS>;
template class KernelArray<TriangularKernel_nxRWG_SS>;

// Quadrilateral instantiation for Planar approximation
template class KernelArray<QuadrilateralKernel_PlanarScalar>;
template class KernelArray<QuadrilateralKernel_PlanarVectorWS>;
template class KernelArray<QuadrilateralKernel_PlanarVectorSS>;

template class KernelArray<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class KernelArray<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class KernelArray<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class KernelArray<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class KernelArray<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class KernelArray<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class KernelArray<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class KernelArray<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class KernelArray<QuadrilateralKernel_CurvilinearScalar>;
template class KernelArray<QuadrilateralKernel_CurvilinearVectorWS>;
template class KernelArray<QuadrilateralKernel_CurvilinearVectorSS>;

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
KernelArray_TR_VA_SS<ParticularKernel>::KernelArray_TR_VA_SS():
KernelArray<ParticularKernel>() {

    // Exclude zero index
    vector<size_t> exclude_idxs = {0};
    this->index_reduction_(exclude_idxs);
}

//virtual
template <typename ParticularKernel>
KernelArray_TR_VA_SS<ParticularKernel>::~KernelArray_TR_VA_SS() {

}

// Triangular elements
template class  KernelArray_TR_VA_SS<TriangularKernel_RWG_SS>;

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
KernelArray_TR_EA_SS<ParticularKernel>::KernelArray_TR_EA_SS():
KernelArray<ParticularKernel>() {

    // Exclude zero index
    vector<size_t> exclude_idxs = {1,3};
    this->index_reduction_(exclude_idxs);
}

//virtual
template <typename ParticularKernel>
KernelArray_TR_EA_SS<ParticularKernel>::~KernelArray_TR_EA_SS() {

}

// Triangular elements
template class KernelArray_TR_EA_SS<TriangularKernel_RWG_SS>;

///////////////////////////////////////////////////////////////////////////////

} // End of the Directfn  namespace

// End of the file
