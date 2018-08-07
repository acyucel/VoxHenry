#include <iostream>
#include "directfn_greenfunc_voxhenry.h"
#include "directfn_common.h"

using  std::cout;
using  std::endl;

namespace Directfn {

//////////////////////////////////////////////////////////////
// VolumeKernel 2
//////////////////////////////////////////////////////////////

VolumeKernel2_KernelType1_GreenFunc::VolumeKernel2_KernelType1_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel2_KernelType1_GreenFunc::~VolumeKernel2_KernelType1_GreenFunc() {
}

//virtual
dcomplex VolumeKernel2_KernelType1_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex K = R / double(2.0) / double(4.0) / M_PI;

    return K;
}

//////////////////////////////////////////////////////////////

VolumeKernel2_KernelType2_GreenFunc::VolumeKernel2_KernelType2_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel2_KernelType2_GreenFunc::~VolumeKernel2_KernelType2_GreenFunc() {
}

//virtual
dcomplex VolumeKernel2_KernelType2_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex K = R / double(8.0) / double(4.0) / M_PI;

    return K;
}

//////////////////////////////////////////////////////////////

VolumeKernel2_KernelType3_GreenFunc::VolumeKernel2_KernelType3_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel2_KernelType3_GreenFunc::~VolumeKernel2_KernelType3_GreenFunc() {
}

//virtual
dcomplex VolumeKernel2_KernelType3_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex K = R / double(8.0) / double(4.0) / M_PI;

    return K;
}

//////////////////////////////////////////////////////////////

VolumeKernel2_KernelType4_GreenFunc::VolumeKernel2_KernelType4_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel2_KernelType4_GreenFunc::~VolumeKernel2_KernelType4_GreenFunc() {
}

//virtual
dcomplex VolumeKernel2_KernelType4_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex K = pow(R,double(3.0)) / double(24.0) / double(4.0) / M_PI;

    return K;
}

//////////////////////////////////////////////////////////////
// VolumeKernel 3
//////////////////////////////////////////////////////////////

VolumeKernel3_KernelType1_GreenFunc::VolumeKernel3_KernelType1_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel3_KernelType1_GreenFunc::~VolumeKernel3_KernelType1_GreenFunc() {

}

//virtual
dcomplex VolumeKernel3_KernelType1_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex K = pow(R,double(3.0)) / double(12.0) * (-k0wn_*k0wn_) / double(8.0) / M_PI;

    return K;
}

//////////////////////////////////////////////////////////////

VolumeKernel3_KernelType2_GreenFunc::VolumeKernel3_KernelType2_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel3_KernelType2_GreenFunc::~VolumeKernel3_KernelType2_GreenFunc() {

}

//virtual
dcomplex VolumeKernel3_KernelType2_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex K = pow(R,double(3.0)) * (-k0wn_*k0wn_) / double(72.0) / double(8.0) / M_PI;

    return K;
}

//////////////////////////////////////////////////////////////

VolumeKernel3_KernelType3_GreenFunc::VolumeKernel3_KernelType3_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel3_KernelType3_GreenFunc::~VolumeKernel3_KernelType3_GreenFunc() {

}

//virtual
dcomplex VolumeKernel3_KernelType3_GreenFunc::genuine_value_(const double R) const noexcept {

	const dcomplex K = pow(R, double(3.0)) * (-k0wn_*k0wn_) / double(72.0) / double(8.0) / M_PI;

    return K;
}

//////////////////////////////////////////////////////////////

VolumeKernel3_KernelType4_GreenFunc::VolumeKernel3_KernelType4_GreenFunc():
AbstractGreenFunc() {
}

//virtual
VolumeKernel3_KernelType4_GreenFunc::~VolumeKernel3_KernelType4_GreenFunc() {

}

//virtual
dcomplex VolumeKernel3_KernelType4_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex K = pow(R,double(5.0)) / double(360.0) * (-k0wn_*k0wn_) / double(8.0) / M_PI;

    return K;
}


}   // End of the namespace Directfn


// End of the file

