#include <iostream>
#include <iomanip>
#include "directfn_kernel_base.h"
#include "directfn_greenfunc.h"

using  std::cout;
using  std::cerr;
using  std::endl;
using  std::setprecision;
using  std::setw;
using  std::to_string;

namespace Directfn {

/////////////////////////////////////////////////////////////////////////////////

AbstractKernel::AbstractKernel():
up_green_func_(nullptr),
rp_crnt_{0.0, 0.0, 0.0},
rq_crnt_{0.0, 0.0, 0.0},
Rpq_{0.0, 0.0, 0.0} {

}

AbstractKernel::~AbstractKernel() {
}

void AbstractKernel::set_wavenumber(const double k0_inp) noexcept {

    if (nullptr == up_green_func_.get()) {
        cout << "AbstractKernel::set_wavenumber Error! " << endl;
        cout << "up_green_func_ is not initialized!... " << endl << endl;
        return;
    }
    up_green_func_->set_wavenumber(k0_inp);
}

void AbstractKernel::precompute_rp_rq_data() noexcept {

    // Rpq_ is used in Green func and in SS cross producrs, so it is not local and prechashed
    for (size_t i = 0; i < 3; ++i) {
        Rpq_[i] = rp_crnt_[i] - rq_crnt_[i];
    }
    up_green_func_->precompute(Rpq_);
    precompute_rp_rq_dependent_data_();
}

dcomplex AbstractKernel::value(const size_t index) const noexcept {
    return specific_value_(index);
}

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel::QuadrilateralKernel():
up_(0.0),
vp_(0.0),
uq_(0.0),
vq_(0.0),
ru_p_{0.0, 0.0, 0.0},
rv_p_{0.0, 0.0, 0.0},
ru_q_{0.0, 0.0, 0.0},
rv_q_{0.0, 0.0, 0.0} {

}

//virtual
QuadrilateralKernel::~QuadrilateralKernel() {
}

double QuadrilateralKernel::calc_Jacobian_() const noexcept {

    double buff[3] = {0.0, 0.0, 0.0};

    vector_cross(ru_p_, rv_p_, buff);
    const double Jp = sqrt(vector_dot(buff, buff));

    vector_cross(ru_q_, rv_q_, buff);
    const double Jq = sqrt(vector_dot(buff, buff));

    return Jp * Jq;
}

///////////////////////////////////////////////////////////////////////////////

} // End of the Directfn  namespace

// End of the file
