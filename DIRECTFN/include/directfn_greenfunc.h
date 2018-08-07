#ifndef _DIRECTFN_GREENFUNC_H_
#define _DIRECTFN_GREENFUNC_H_

#include "directfn_defs.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

class AbstractGreenFunc {
public:

    AbstractGreenFunc();
    ~AbstractGreenFunc();

    AbstractGreenFunc(const AbstractGreenFunc & ) = delete;
    AbstractGreenFunc(AbstractGreenFunc && ) = delete;
    AbstractGreenFunc & operator = (const AbstractGreenFunc & ) = delete;
    AbstractGreenFunc & operator = (AbstractGreenFunc && ) = delete;

    /*! The wave number is used in the exp(i k R)  */
    void set_wavenumber(const double k0_inp) noexcept;

    /*! Perfoms computation of the Green function value and saves it.
     *  In fact, it calls the genuine_value_() routine. */
    void  precompute(const double Rpq[3]) noexcept;

    /*! Return precomputed value. */
    dcomplex  value() const noexcept;

    virtual void debug_print() const noexcept;

protected:
    /*! The wave number is setup once and for all (but can be changed if needed) */
    double k0wn_;

private:
    /*! Saves the result here */
    dcomplex precomputed_value_;

    /*! Calls the actual value */
    virtual dcomplex genuine_value_(const double R) const noexcept = 0;
};

///////////////////////////////////////////////////////////////////////////////

class HelmgolzGreenFunc final : public AbstractGreenFunc {
public:
    HelmgolzGreenFunc();
    ~HelmgolzGreenFunc();

private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class GradHelmgolzGreenFunc final : public AbstractGreenFunc {
public:
    GradHelmgolzGreenFunc();
    ~GradHelmgolzGreenFunc();

private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class G_minus_G0_GreenFunc final : public AbstractGreenFunc {
public:
    G_minus_G0_GreenFunc();
    virtual ~G_minus_G0_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class Grad_G_minus_G0_GreenFunc final : public AbstractGreenFunc {
public:
    Grad_G_minus_G0_GreenFunc();
    virtual ~Grad_G_minus_G0_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

class Grad_G_minus_half_G0_GreenFunc final : public AbstractGreenFunc {
public:
    Grad_G_minus_half_G0_GreenFunc();
    virtual ~Grad_G_minus_half_G0_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace  Directfn

#endif   //  _DIRECTFN_GREENFUNC_H_

// End of the file

