#ifndef _DIRECTFN_GREENFUNC_VOXHENRY_H_
#define _DIRECTFN_GREENFUNC_VOXHENRY_H_

#include "directfn_defs.h"
#include "directfn_greenfunc.h"

namespace Directfn {


///////////////////////////////////////////////////////////////////////////////

class VolumeKernel2_KernelType1_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel2_KernelType1_GreenFunc();
    virtual ~VolumeKernel2_KernelType1_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

class VolumeKernel2_KernelType2_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel2_KernelType2_GreenFunc();
    virtual ~VolumeKernel2_KernelType2_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

class VolumeKernel2_KernelType3_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel2_KernelType3_GreenFunc();
    virtual ~VolumeKernel2_KernelType3_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

class VolumeKernel2_KernelType4_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel2_KernelType4_GreenFunc();
    virtual ~VolumeKernel2_KernelType4_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

class VolumeKernel3_KernelType2_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel3_KernelType2_GreenFunc();
    virtual ~VolumeKernel3_KernelType2_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

class VolumeKernel3_KernelType3_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel3_KernelType3_GreenFunc();
    virtual ~VolumeKernel3_KernelType3_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

class VolumeKernel3_KernelType1_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel3_KernelType1_GreenFunc();
    virtual ~VolumeKernel3_KernelType1_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};

class VolumeKernel3_KernelType4_GreenFunc final : public AbstractGreenFunc {
public:
    VolumeKernel3_KernelType4_GreenFunc();
    virtual ~VolumeKernel3_KernelType4_GreenFunc();
private:
    virtual dcomplex genuine_value_(const double R) const noexcept;
};


///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace  Directfn

#endif   //  _DIRECTFN_GREENFUNC_VOXHENRY_H_

// End of the file

