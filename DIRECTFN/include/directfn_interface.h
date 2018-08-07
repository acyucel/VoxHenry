#ifndef _DIRECTFN_INTERFACE_H_
#define _DIRECTFN_INTERFACE_H_

#include <memory>
#include <vector>
#include <string>
#include "directfn_kernel_arrsum.h"
#include "directfn_contour.h"
#include "directfn_defs.h"

using std::string;
using std::unique_ptr;
using std::shared_ptr;
using std::size_t;
using std::vector;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

class SingularContour3xn;

/*! \class DirectfnInterface  ....
 *
 *  Interface for directfn
 */
template <typename ParticularKernel>
class DirectfnInterface {
public:
    /*! Default constructor: all containers are of zero size,
     *  all pointers are nullptr initialized. */
    DirectfnInterface();

    /*! Default destructor */
    virtual ~DirectfnInterface();

    /*! The internal data are not supposed to be copied neither in the
     *  copy constructor nor in the assignment operator. */
    DirectfnInterface(const DirectfnInterface & ) = delete;
    DirectfnInterface(DirectfnInterface && ) = delete;
    DirectfnInterface & operator = (const DirectfnInterface & ) = delete;
    DirectfnInterface & operator = (DirectfnInterface && ) = delete;

    /*! Wave number setup for kernel exponent (k * R * Imaginary_unity) */
    virtual void set_wavenumber(const double k0_inp) noexcept;

    /*! The orders of Gaussian quadratures in each dimension in the transformed (final) quadrangle(s). */
    bool set_Gaussian_orders_4(const size_t N1, const size_t N2, const size_t N3, const size_t N4) noexcept;
    bool set_Gaussian_orders_4(const size_t Nx[4]) noexcept;
    bool set_Gaussian_orders_4(const size_t N_the_same_1234) noexcept;

    void set(const SingularContour3xn &  contour_xpts) noexcept;

    size_t kernel_size() const noexcept;

    size_t  N1() const noexcept;
    size_t  N2() const noexcept;
    size_t  N3() const noexcept;
    size_t  N4() const noexcept;

    /*! Peforms summation for all integrals between the kernel and all testing and basis functions. */
    void calc_I_surface_surface();
    void calc_Iss();

    /*! Ref to array of kernel-I_surf_surf */
    const dcomplex * Iss()     const noexcept;
    const dcomplex   Iss_tot() const noexcept;
    const dcomplex   Iss_arr(const size_t k) const noexcept;

    void copy_Iss_array_values_to(dcomplex * const out_array_allocated) const noexcept;

    /*! \return kernel pointer to setup parametric kernel properly */
    ParticularKernel * kernel_ptr() noexcept;

    /*! The constant string name of the kernel needed usually for debug */
    virtual string name() const noexcept = 0;

    virtual void debug_print() const noexcept;

protected:
    /*! Kernels pointers will be allocated in the specialized constructors */
    unique_ptr<ParticularKernel> up_kernel_;

    /*! Summarize integral values over rho_4, lam_3, eta_2 in kernel_sz loop */
    unique_ptr<KernelArrayInterface<ParticularKernel>>  up_kerSummator_;

    /*! Pointers to arrays of points and weights for Gaussian quadratures. */
    unique_ptr<double []>  up_w1_;
    unique_ptr<double []>  up_z1_;
    unique_ptr<double []>  up_w2_;
    unique_ptr<double []>  up_z2_;
    unique_ptr<double []>  up_w3_;
    unique_ptr<double []>  up_z3_;
    unique_ptr<double []>  up_w4_;
    unique_ptr<double []>  up_z4_;

    /*! Surface-surface values */
    unique_ptr<dcomplex []>  Iss_;

    /*! Iss sum over all kernels */
    dcomplex Iss_tot_;

private:
    /*! Orders of Gaussian quadratures */
    size_t  N1_, N2_, N3_, N4_;

    /*! Main integration routine redefined in inherited classes */
    virtual void do_I_surface_surface_() = 0;

    /*! Auxiliarily routine for Gauss-Legandre nodes setup */
    bool set_zw_N_(const size_t Nx, unique_ptr<double []> & p_z1, unique_ptr<double []> & p_w1);
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Dierectfn

#endif  // _DIRECTFN_INTERFACE_H_

// End of the file




