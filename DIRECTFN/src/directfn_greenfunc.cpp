#include <iostream>
#include "directfn_greenfunc.h"
#include "directfn_common.h"

using  std::cout;
using  std::endl;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

AbstractGreenFunc::AbstractGreenFunc():
k0wn_(0.0),
precomputed_value_(0.0, 0.0) {

}

AbstractGreenFunc::~AbstractGreenFunc() {
}

void AbstractGreenFunc::set_wavenumber(const double k0_inp) noexcept {
    k0wn_ = k0_inp;
}

void AbstractGreenFunc::precompute(const double Rpq[3]) noexcept {

    const double tR = sqrt(vector_dot(Rpq, Rpq));
    precomputed_value_ = genuine_value_(tR);
}

dcomplex AbstractGreenFunc::value() const noexcept {
    return precomputed_value_;
}

//virtual
void AbstractGreenFunc::debug_print() const noexcept {
    cout << "AbstractGreenFunc: k0 = " << k0wn_ << endl;
}

///////////////////////////////////////////////////////////////////////////////

HelmgolzGreenFunc::HelmgolzGreenFunc():
AbstractGreenFunc() {
}

HelmgolzGreenFunc::~HelmgolzGreenFunc() {
}

//virtual
dcomplex HelmgolzGreenFunc::genuine_value_(const double R) const noexcept {
    return exp(-Iunit * k0wn_ * R) / (double(4.0) * M_PI * R);
}

///////////////////////////////////////////////////////////////////////////////

GradHelmgolzGreenFunc::GradHelmgolzGreenFunc():
AbstractGreenFunc() {
}

GradHelmgolzGreenFunc::~GradHelmgolzGreenFunc() {
}

//virtual
dcomplex GradHelmgolzGreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex Green = exp(-Iunit * k0wn_ * R) / (4.0 * M_PI * R);
    return  -Green * (Iunit * k0wn_ * R + 1.0) / (R * R);
}

///////////////////////////////////////////////////////////////////////////////

G_minus_G0_GreenFunc::G_minus_G0_GreenFunc():
AbstractGreenFunc() {
}

//virtual
G_minus_G0_GreenFunc::~G_minus_G0_GreenFunc() {
}

//virtual
dcomplex G_minus_G0_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex Green_0 =  1.0 / (4.0 * M_PI * R);
    const dcomplex Green   = exp(-Iunit * k0wn_ * R) * Green_0;

    return (Green_0 - Green) / (k0wn_ * k0wn_);
}

///////////////////////////////////////////////////////////////////////////////

Grad_G_minus_G0_GreenFunc::Grad_G_minus_G0_GreenFunc():
AbstractGreenFunc() {
}

//virtual
Grad_G_minus_G0_GreenFunc::~Grad_G_minus_G0_GreenFunc() {

}

//virtual
dcomplex Grad_G_minus_G0_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex G0 = 1.0 / (4.0 * M_PI * R);
    const dcomplex G  = exp(-Iunit * k0wn_ * R) * G0;

    const dcomplex GRfactor = (-Iunit * k0wn_ * G * R - G + G0) / (R * R * k0wn_ * k0wn_);

    return GRfactor;
}

///////////////////////////////////////////////////////////////////////////////

Grad_G_minus_half_G0_GreenFunc::Grad_G_minus_half_G0_GreenFunc():
AbstractGreenFunc() {
}

Grad_G_minus_half_G0_GreenFunc::~Grad_G_minus_half_G0_GreenFunc() {
}

//virtual
dcomplex Grad_G_minus_half_G0_GreenFunc::genuine_value_(const double R) const noexcept {

    const dcomplex  G0 = 1.0 / (4.0 * M_PI * R);
    const dcomplex  G  = exp(-Iunit * k0wn_ * R) * G0;

    const dcomplex  GRfactor = (-Iunit * k0wn_ * G * R - G + G0) / (R * R * k0wn_ * k0wn_);

    return GRfactor +  0.5 * G0;
}

///////////////////////////////////////////////////////////////////////////////

}   // End of the namespace Directfn

// End of the file

