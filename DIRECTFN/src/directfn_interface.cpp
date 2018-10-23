#include <iostream>
#include <iomanip>
#include "directfn_interface.h"
#include "directfn_kernel_tri.h"
#include "directfn_kernel_quad_scal.h"
#include "directfn_kernel_quad_vect.h"
#include "directfn_kernel_quad_voxhenry.h"

using std::cout;
using std::endl;
using std::setw;


namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
DirectfnInterface<ParticularKernel>::DirectfnInterface():
up_kernel_(nullptr),
up_kerSummator_(nullptr),
up_w1_(nullptr),
up_z1_(nullptr),
up_w2_(nullptr),
up_z2_(nullptr),
up_w3_(nullptr),
up_z3_(nullptr),
up_w4_(nullptr),
up_z4_(nullptr),
Iss_(nullptr),
Iss_tot_(0.0, 0.0),
N1_(0),
N2_(0),
N3_(0),
N4_(0) {

    up_kernel_.reset(new ParticularKernel());
    up_kerSummator_.reset(new KernelArray<ParticularKernel>());
    up_kerSummator_->setup(up_kernel_.get());
}

//virtual
template <typename ParticularKernel>
DirectfnInterface<ParticularKernel>::~DirectfnInterface() {
}

template <typename ParticularKernel>
void DirectfnInterface<ParticularKernel>::set_wavenumber(const double k0_inp) noexcept {
    up_kernel_->set_wavenumber(k0_inp);
}

template <typename ParticularKernel>
bool DirectfnInterface<ParticularKernel>::set_Gaussian_orders_4(const size_t N1_in, const size_t N2_in,
                                                                const size_t N3_in, const size_t N4_in) noexcept {
    N1_ = N1_in;
    N2_ = N2_in;
    N3_ = N3_in;
    N4_ = N4_in;

    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N1_, up_z1_, up_w1_)) {return false;}
    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N2_, up_z2_, up_w2_)) {return false;}
    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N3_, up_z3_, up_w3_)) {return false;}
    if (!DirectfnInterface<ParticularKernel>::set_zw_N_(N4_, up_z4_, up_w4_)) {return false;}

    return true;
}

template <typename ParticularKernel>
bool DirectfnInterface<ParticularKernel>::set_Gaussian_orders_4(const size_t Nx[4]) noexcept {
    return DirectfnInterface::set_Gaussian_orders_4(Nx[0], Nx[1], Nx[2], Nx[3]);
}

template <typename ParticularKernel>
bool DirectfnInterface<ParticularKernel>::set_Gaussian_orders_4(const size_t Np) noexcept {
    return this->set_Gaussian_orders_4(Np, Np, Np, Np);
}

template <typename ParticularKernel>
void DirectfnInterface<ParticularKernel>::set(const SingularContour3xn &  contour_xpts) noexcept {
    up_kernel_->set(contour_xpts);
}

template <typename ParticularKernel>
size_t DirectfnInterface<ParticularKernel>::kernel_size() const noexcept {
    return up_kernel_->size();
}

template <typename ParticularKernel>
size_t  DirectfnInterface<ParticularKernel>::N1() const noexcept {
    return N1_;
}

template <typename ParticularKernel>
size_t  DirectfnInterface<ParticularKernel>::N2() const noexcept {
    return N2_;
}

template <typename ParticularKernel>
size_t  DirectfnInterface<ParticularKernel>::N3() const noexcept {
    return N3_;
}

template <typename ParticularKernel>
size_t  DirectfnInterface<ParticularKernel>::N4() const noexcept {
    return N4_;
}

template <typename ParticularKernel>
void DirectfnInterface<ParticularKernel>::calc_I_surface_surface() {
    return this->do_I_surface_surface_();
}

template <typename ParticularKernel>
void DirectfnInterface<ParticularKernel>::calc_Iss() {
    return this->do_I_surface_surface_();
}

template <typename ParticularKernel>
const dcomplex * DirectfnInterface<ParticularKernel>::Iss() const noexcept {
    return Iss_.get();
}

template <typename ParticularKernel>
const dcomplex DirectfnInterface<ParticularKernel>::Iss_tot() const noexcept {
    return Iss_tot_;
}

template <typename ParticularKernel>
const dcomplex DirectfnInterface<ParticularKernel>::Iss_arr(const size_t k) const noexcept {
    return Iss_[k];
}

template <typename ParticularKernel>
void DirectfnInterface<ParticularKernel>::copy_Iss_array_values_to(dcomplex * const out_array_to_be_setup) const noexcept {

    for (size_t k = 0; k < this->kernel_size(); ++k) {
        out_array_to_be_setup[k] = Iss_[k];
    }
}

template <typename ParticularKernel>
ParticularKernel * DirectfnInterface<ParticularKernel>::kernel_ptr() noexcept {
    return up_kernel_.get();
}

//virtual
template <typename ParticularKernel>
void DirectfnInterface<ParticularKernel>::debug_print() const noexcept {
    cout << "N1 = " << N1_ << "   N2 = " << N2_ << "   N3_ = " << N3_ << "   N4 = " << N4_ << endl;
}

template <typename ParticularKernel>
bool DirectfnInterface<ParticularKernel>::set_zw_N_(const size_t Nn,
                                                    unique_ptr<double []> & up_zn,
                                                    unique_ptr<double []> & up_wn) {
    up_zn.reset(new double[Nn]);
    up_wn.reset(new double[Nn]);

    gl_xw_1d(int(Nn), up_zn.get(), up_wn.get());
    return true;
}


// Instantiation of Triangular elements
template class DirectfnInterface<TriangularKernel_Constant_ST>;
template class DirectfnInterface<TriangularKernel_Constant_EA>;
template class DirectfnInterface<TriangularKernel_Constant_VA>;
template class DirectfnInterface<TriangularKernel_RWG_WS>;
template class DirectfnInterface<TriangularKernel_RWG_SS>;
template class DirectfnInterface<TriangularKernel_nxRWG_SS>;

// Instantiation of Quadrilateral Constant kernels
template class DirectfnInterface<QuadrilateralKernel_PlanarScalar>;
template class DirectfnInterface<QuadrilateralKernel_PlanarVectorWS>;
template class DirectfnInterface<QuadrilateralKernel_PlanarVectorSS>;

template class DirectfnInterface<QuadKer_PlanVH_VolKer2_KerTyp1>;
template class DirectfnInterface<QuadKer_PlanVH_VolKer2_KerTyp2>;
template class DirectfnInterface<QuadKer_PlanVH_VolKer2_KerTyp3>;
template class DirectfnInterface<QuadKer_PlanVH_VolKer2_KerTyp4>;
template class DirectfnInterface<QuadKer_PlanVH_VolKer3_KerTyp1>;
template class DirectfnInterface<QuadKer_PlanVH_VolKer3_KerTyp2>;
template class DirectfnInterface<QuadKer_PlanVH_VolKer3_KerTyp3>;
template class DirectfnInterface<QuadKer_PlanVH_VolKer3_KerTyp4>;

template class DirectfnInterface<QuadrilateralKernel_CurvilinearScalar>;
template class DirectfnInterface<QuadrilateralKernel_CurvilinearVectorWS>;
template class DirectfnInterface<QuadrilateralKernel_CurvilinearVectorSS>;

///////////////////////////////////////////////////////////////////////////////

}   // namespace Directfn

// End of the file


