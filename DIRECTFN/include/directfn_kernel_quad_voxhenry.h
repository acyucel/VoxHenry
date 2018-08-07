#ifndef _DIRECTFN_KERNEL_QUAD_VOXHENRY_H_
#define _DIRECTFN_KERNEL_QUAD_VOXHENRY_H_

#include "directfn_kernel_base.h"
#include "directfn_kernel_quad_geom.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_VoxHenry : public virtual QuadrilateralKernel {
public:
    QuadrilateralKernel_VoxHenry();
    virtual ~QuadrilateralKernel_VoxHenry();

    QuadrilateralKernel_VoxHenry(const QuadrilateralKernel_VoxHenry & ) = delete;
    QuadrilateralKernel_VoxHenry(QuadrilateralKernel_VoxHenry && ) = delete;
    QuadrilateralKernel_VoxHenry & operator = (const QuadrilateralKernel_VoxHenry & ) = delete;
    QuadrilateralKernel_VoxHenry & operator = (QuadrilateralKernel_VoxHenry && ) = delete;
   
    void set_lp(const int lp) noexcept;
    
protected:
    /*! lp index */
    int lp_;
};

class QuadrilateralPlanarKernel_VoxHenry : public virtual QuadrilateralPlanarKernel {
public:
    QuadrilateralPlanarKernel_VoxHenry();
    virtual ~QuadrilateralPlanarKernel_VoxHenry();

    QuadrilateralPlanarKernel_VoxHenry(const QuadrilateralPlanarKernel_VoxHenry & ) = delete;
    QuadrilateralPlanarKernel_VoxHenry(QuadrilateralPlanarKernel_VoxHenry && ) = delete;
    QuadrilateralPlanarKernel_VoxHenry & operator = (const QuadrilateralPlanarKernel_VoxHenry & ) = delete;
    QuadrilateralPlanarKernel_VoxHenry & operator = (QuadrilateralPlanarKernel_VoxHenry & ) = delete;

    void set_nq_np(const double nq[3], const double np[3]) noexcept;
    
protected:
    /*! normals to the quadrilaterals */
    double np_[3], nq_[3];
    
    /*! side size */
    double dx_;
    
    /*! np_, nq_ i.e. normals to the quadrilaterals */
    virtual void   precompute_vertex_dependent_data_() noexcept;
};


///////////////////////////////////////////////////////////////////////////////

class QuadrilateralKernel_PlanarVoxHenry : public QuadrilateralPlanarKernel_VoxHenry,
                                               public QuadrilateralKernel_VoxHenry {
public:
    QuadrilateralKernel_PlanarVoxHenry();
    virtual ~QuadrilateralKernel_PlanarVoxHenry();

    QuadrilateralKernel_PlanarVoxHenry(const QuadrilateralKernel_PlanarVoxHenry & ) = delete;
    QuadrilateralKernel_PlanarVoxHenry(QuadrilateralKernel_PlanarVoxHenry && ) = delete;

    QuadrilateralKernel_PlanarVoxHenry & operator = (const QuadrilateralKernel_PlanarVoxHenry & ) = delete;
    QuadrilateralKernel_PlanarVoxHenry & operator = (QuadrilateralKernel_PlanarVoxHenry && ) = delete;

    virtual size_t size() const noexcept;
    static constexpr size_t constexpr_size() noexcept {return 1;}

    /*! Is used after N1...N4 loop computations and needed if the Jacobian is constant.
     *  Returns 1 because it varies inside the kernel. */
    virtual double precomputed_jacobian() const noexcept;
    
protected:

    void calculate_Np_() noexcept;
    void calculate_Nq_() noexcept;
    
    // Nodal shape functions
	double Np_;  
	double Nq_;
	
	// dot prodct part of the kernel
	double F_;
};



class QuadKer_PlanVH_VolKer2_KerTyp1 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer2_KerTyp1();
    virtual ~QuadKer_PlanVH_VolKer2_KerTyp1();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

class QuadKer_PlanVH_VolKer2_KerTyp2 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer2_KerTyp2();
    virtual ~QuadKer_PlanVH_VolKer2_KerTyp2();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

class QuadKer_PlanVH_VolKer2_KerTyp3 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer2_KerTyp3();
    virtual ~QuadKer_PlanVH_VolKer2_KerTyp3();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

class QuadKer_PlanVH_VolKer2_KerTyp4 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer2_KerTyp4();
    virtual ~QuadKer_PlanVH_VolKer2_KerTyp4();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

class QuadKer_PlanVH_VolKer3_KerTyp1 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer3_KerTyp1();
    virtual ~QuadKer_PlanVH_VolKer3_KerTyp1();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

class QuadKer_PlanVH_VolKer3_KerTyp2 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer3_KerTyp2();
    virtual ~QuadKer_PlanVH_VolKer3_KerTyp2();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

class QuadKer_PlanVH_VolKer3_KerTyp3 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer3_KerTyp3();
    virtual ~QuadKer_PlanVH_VolKer3_KerTyp3();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

class QuadKer_PlanVH_VolKer3_KerTyp4 final : public QuadrilateralKernel_PlanarVoxHenry {
public:
    QuadKer_PlanVH_VolKer3_KerTyp4();
    virtual ~QuadKer_PlanVH_VolKer3_KerTyp4();

protected:
    virtual void precompute_rp_rq_dependent_data_()  noexcept;

private:
    virtual dcomplex specific_value_(const size_t ) const noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  //  _DIRECTFN_KERNEL_QUAD_VOXHENRY_H_

// End of the file
