#include <iostream>

#include "directfn_kernel_quad_voxhenry.h"
#include "directfn_greenfunc_voxhenry.h"

using std::cout;
using std::endl;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_VoxHenry::QuadrilateralKernel_VoxHenry():
QuadrilateralKernel() {
}

//virtual
QuadrilateralKernel_VoxHenry::~QuadrilateralKernel_VoxHenry() {
}

void QuadrilateralKernel_VoxHenry::set_lp(const int lp) noexcept {
    lp_ = lp;
}


///////////////////////////////////////////////////////////////////////////////

QuadrilateralPlanarKernel_VoxHenry::QuadrilateralPlanarKernel_VoxHenry():
QuadrilateralPlanarKernel() {

}

//virtual
QuadrilateralPlanarKernel_VoxHenry::~QuadrilateralPlanarKernel_VoxHenry() {

}

//virtual
void QuadrilateralPlanarKernel_VoxHenry::precompute_vertex_dependent_data_() noexcept {
   
    QuadrilateralPlanarKernel::precompute_vertex_dependent_data_();

/*  Commented out, as normals to p and q must be set from the external,
    and not be inferred by the vertex order, as in the case of 6 points
    this cannot be done.

	// Compute normal to the quadrilaterals, using Newell's method (see e.g. graphics gems III, V.5)
	// Normal direction is so that N x rp1_-rp4_
	// is directed inside the polygon

    double t_norm;    

	// calculate components of normal vector, using the
	// fact that the area of the trapezoids projected on the
	// xy, xz, yz planes is proportional to the z, y, x
	// components of the normal

    // normal to quadrilateral p
    //
    
    np_[0] = (rp1_[1] - rp2_[1]) * (rp1_[2] + rp2_[2]);
    np_[1] = (rp1_[2] - rp2_[2]) * (rp1_[0] + rp2_[0]);
    np_[2] = (rp1_[0] - rp2_[0]) * (rp1_[1] + rp2_[1]);
    
    np_[0] += (rp2_[1] - rp3_[1]) * (rp2_[2] + rp3_[2]);
    np_[1] += (rp2_[2] - rp3_[2]) * (rp2_[0] + rp3_[0]);
    np_[2] += (rp2_[0] - rp3_[0]) * (rp2_[1] + rp3_[1]);

    np_[0] += (rp3_[1] - rp4_[1]) * (rp3_[2] + rp4_[2]);
    np_[1] += (rp3_[2] - rp4_[2]) * (rp3_[0] + rp4_[0]);
    np_[2] += (rp3_[0] - rp4_[0]) * (rp3_[1] + rp4_[1]);
   
    np_[0] += (rp4_[1] - rp1_[1]) * (rp4_[2] + rp1_[2]);
    np_[1] += (rp4_[2] - rp1_[2]) * (rp4_[0] + rp1_[0]);
    np_[2] += (rp4_[0] - rp1_[0]) * (rp4_[1] + rp1_[1]);

    t_norm = sqrt(vector_dot(np_, np_));
    
    np_[0] /= t_norm;
    np_[1] /= t_norm;
    np_[2] /= t_norm;
   
    // normal to quadrilateral q
    //
    
    nq_[0] = (rq1_[1] - rq2_[1]) * (rq1_[2] + rq2_[2]);
    nq_[1] = (rq1_[2] - rq2_[2]) * (rq1_[0] + rq2_[0]);
    nq_[2] = (rq1_[0] - rq2_[0]) * (rq1_[1] + rq2_[1]);
    
    nq_[0] += (rq2_[1] - rq3_[1]) * (rq2_[2] + rq3_[2]);
    nq_[1] += (rq2_[2] - rq3_[2]) * (rq2_[0] + rq3_[0]);
    nq_[2] += (rq2_[0] - rq3_[0]) * (rq2_[1] + rq3_[1]);

    nq_[0] += (rq3_[1] - rq4_[1]) * (rq3_[2] + rq4_[2]);
    nq_[1] += (rq3_[2] - rq4_[2]) * (rq3_[0] + rq4_[0]);
    nq_[2] += (rq3_[0] - rq4_[0]) * (rq3_[1] + rq4_[1]);
   
    nq_[0] += (rq4_[1] - rq1_[1]) * (rq4_[2] + rq1_[2]);
    nq_[1] += (rq4_[2] - rq1_[2]) * (rq4_[0] + rq1_[0]);
    nq_[2] += (rq4_[0] - rq1_[0]) * (rq4_[1] + rq1_[1]);

    t_norm = sqrt(vector_dot(nq_, nq_));
    
    nq_[0] /= t_norm;
    nq_[1] /= t_norm;
    nq_[2] /= t_norm;
*/   

    // size of the side
    double lp[3];
    lp[0] = (rp2_[0] - rp1_[0]);
    lp[1] = (rp2_[1] - rp1_[1]);
    lp[2] = (rp2_[2] - rp1_[2]);
    // warning: should check that all sides of both q and p have the same dimensions!
    dx_ = sqrt(vector_dot(lp, lp));
    
    // debug
    
    //cout << "lp " << lp[0] << ' ' << lp[1] << ' ' << lp[2] << endl;
    //cout << "dx_ " << dx_  << endl;
    //cout << "np " << np_[0] << ' ' << np_[1] << ' ' << np_[2] << endl;
    //cout << "nq " << nq_[0] << ' ' << nq_[1] << ' ' << nq_[2] << endl;
    //cout << "cp " << rp_pppp_[0] << ' ' << rp_pppp_[1] << ' ' << rp_pppp_[2] << endl;
    //cout << "cq " << rq_pppp_[0] << ' ' << rq_pppp_[1] << ' ' << rq_pppp_[2] << endl;

}

void QuadrilateralPlanarKernel_VoxHenry::set_nq_np(const double nq[3], const double np[3]) noexcept {
    int i;
    
    for(i=0; i<3; i++) {
        nq_[i] = nq[i];
        np_[i] = np[i];
    }
}

///////////////////////////////////////////////////////////////////////////////

QuadrilateralKernel_PlanarVoxHenry::QuadrilateralKernel_PlanarVoxHenry():
QuadrilateralKernel(),
QuadrilateralPlanarKernel_VoxHenry(),
QuadrilateralKernel_VoxHenry() {

}


// virtual
QuadrilateralKernel_PlanarVoxHenry::~QuadrilateralKernel_PlanarVoxHenry() {

}

//virtual
size_t QuadrilateralKernel_PlanarVoxHenry::size() const noexcept {
    return constexpr_size();
}

//virtual
double QuadrilateralKernel_PlanarVoxHenry::precomputed_jacobian() const noexcept {
    // This is not jacobian_ used in the loop where rp_ rq_ changes.
    // It must be used if and only if the J is const for all rp and rq
    
    double jacobian = calc_Jacobian_();
    
    //cout << "jacobian " << jacobian << endl;
    
    return jacobian;
}


void QuadrilateralKernel_PlanarVoxHenry::calculate_Np_() noexcept {

    // rp_pppp_, rq_pppp_ are actually the centers of the quadrilaterals
    // np_, nq_ are the quadrilateral normals
    // rp_crnt_, rq_crnt_ are the coordinates of the numerical integration points
    if (lp_ == 0)
    {
        Np_ = double(1.0);
    }
    else if (lp_ == 1)
    {
        Np_ = double(1.0);
    }
    else if (lp_ == 2)
    {
        Np_ = rp_crnt_[0] - rp_pppp_[0] + double(0.5)*np_[0] * dx_;
    }
    else if (lp_ == 3)
    {
        Np_ = double(1.0);
    }
    else if (lp_ == 4)
    {
        Np_ = rp_crnt_[1] - rp_pppp_[1] + double(0.5)*np_[1] * dx_;
    }
    else if (lp_ == 5)
    {
        Np_ = double(1.0);
    }
    else if (lp_ == 6)
    {	
        Np_ = rp_crnt_[2] - rp_pppp_[2] + double(0.5)*np_[2] * dx_;
    }
    else if (lp_ == 7)
    {
        Np_ = rp_crnt_[0] - rp_pppp_[0] + double(0.5)*np_[0] * dx_;
    }
    else if (lp_ == 8)
    {
        Np_ = rp_crnt_[1] - rp_pppp_[1] + double(0.5)*np_[1] * dx_;
    }
    else if (lp_ == 9)
    {	
        Np_ = rp_crnt_[2] - rp_pppp_[2] + double(0.5)*np_[2] * dx_;
    }
    else
    {
        // this should never happen
        Np_ = double(1.0);
    }
}

void QuadrilateralKernel_PlanarVoxHenry::calculate_Nq_() noexcept {

    // rp_pppp_, rq_pppp_ are actually the centers of the quadrilaterals
    // np_, nq_ are the quadrilateral normals
    // rp_crnt_, rq_crnt_ are the coordinates of the numerical integration points
    if (lp_ == 0)
    {
        Nq_ = double(1.0);
    }
    else if (lp_ == 1)
    {
        Nq_ = rq_crnt_[0] - rq_pppp_[0] + double(0.5)*nq_[0] * dx_;
    }
    else if (lp_ == 2)
    {
        Nq_ = double(1.0);
    }
    else if (lp_ == 3)
    {
        Nq_ = rq_crnt_[1] - rq_pppp_[1] + double(0.5)*nq_[1] * dx_;
    }
    else if (lp_ == 4)
    {
        Nq_ = double(1.0);
    }
    else if (lp_ == 5)
    {
        Nq_ = rq_crnt_[2] - rq_pppp_[2] + double(0.5)*nq_[2] * dx_;
    }
    else if (lp_ == 6)
    {	
        Nq_ = double(1.0);
    }
    else if (lp_ == 7)
    {
        Nq_ = rq_crnt_[0] - rq_pppp_[0] + double(0.5)*nq_[0] * dx_;
    }
    else if (lp_ == 8)
    {
        Nq_ = rq_crnt_[1] - rq_pppp_[1] + double(0.5)*nq_[1] * dx_;
    }
    else if (lp_ == 9)
    {	
        Nq_ = rq_crnt_[2] - rq_pppp_[2] + double(0.5)*nq_[2] * dx_;
    }
    else
    {
        // this should never happen
        Nq_ = double(1.0);
    }
}


///////////////////////////////////////////////////////////////////////////////
// VolumeKernel 2
///////////////////////////////////////////////////////////////////////////////

QuadKer_PlanVH_VolKer2_KerTyp1::QuadKer_PlanVH_VolKer2_KerTyp1() {
    up_green_func_.reset(new VolumeKernel2_KernelType1_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer2_KerTyp1::~QuadKer_PlanVH_VolKer2_KerTyp1() {
}

// virtual
void QuadKer_PlanVH_VolKer2_KerTyp1::precompute_rp_rq_dependent_data_()  noexcept {
    
    calculate_Np_();
    calculate_Nq_();
}


//virtual
dcomplex QuadKer_PlanVH_VolKer2_KerTyp1::specific_value_(const size_t ) const noexcept {
    
    return Np_*Nq_ * up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////

QuadKer_PlanVH_VolKer2_KerTyp2::QuadKer_PlanVH_VolKer2_KerTyp2() {
    up_green_func_.reset(new VolumeKernel2_KernelType2_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer2_KerTyp2::~QuadKer_PlanVH_VolKer2_KerTyp2() {
}

// virtual
void QuadKer_PlanVH_VolKer2_KerTyp2::precompute_rp_rq_dependent_data_()  noexcept {
    
    calculate_Nq_();
       
    F_ = vector_dot(np_, Rpq_);
    
}


//virtual
dcomplex QuadKer_PlanVH_VolKer2_KerTyp2::specific_value_(const size_t ) const noexcept {
    
    return Nq_ * F_ * up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////

QuadKer_PlanVH_VolKer2_KerTyp3::QuadKer_PlanVH_VolKer2_KerTyp3() {
    up_green_func_.reset(new VolumeKernel2_KernelType3_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer2_KerTyp3::~QuadKer_PlanVH_VolKer2_KerTyp3() {
}

// virtual
void QuadKer_PlanVH_VolKer2_KerTyp3::precompute_rp_rq_dependent_data_()  noexcept {
    
    calculate_Np_();
    
    F_ = vector_dot(np_, Rpq_);
}


//virtual
dcomplex QuadKer_PlanVH_VolKer2_KerTyp3::specific_value_(const size_t ) const noexcept {
    
    return Np_ * F_ * up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////


QuadKer_PlanVH_VolKer2_KerTyp4::QuadKer_PlanVH_VolKer2_KerTyp4() {
    up_green_func_.reset(new VolumeKernel2_KernelType4_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer2_KerTyp4::~QuadKer_PlanVH_VolKer2_KerTyp4() {
}

// virtual
void QuadKer_PlanVH_VolKer2_KerTyp4::precompute_rp_rq_dependent_data_()  noexcept {
}

//virtual
dcomplex QuadKer_PlanVH_VolKer2_KerTyp4::specific_value_(const size_t ) const noexcept {
    return up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////
// VolumeKernel 3
///////////////////////////////////////////////////////////////////////////////


QuadKer_PlanVH_VolKer3_KerTyp1::QuadKer_PlanVH_VolKer3_KerTyp1() {
    up_green_func_.reset(new VolumeKernel3_KernelType1_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer3_KerTyp1::~QuadKer_PlanVH_VolKer3_KerTyp1() {
}

// virtual
void QuadKer_PlanVH_VolKer3_KerTyp1::precompute_rp_rq_dependent_data_()  noexcept {
    
    calculate_Np_();
    calculate_Nq_();
}


//virtual
dcomplex QuadKer_PlanVH_VolKer3_KerTyp1::specific_value_(const size_t ) const noexcept {
    
    return Np_*Nq_ * up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////


QuadKer_PlanVH_VolKer3_KerTyp2::QuadKer_PlanVH_VolKer3_KerTyp2() {
    up_green_func_.reset(new VolumeKernel3_KernelType2_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer3_KerTyp2::~QuadKer_PlanVH_VolKer3_KerTyp2() {
}

// virtual
void QuadKer_PlanVH_VolKer3_KerTyp2::precompute_rp_rq_dependent_data_()  noexcept {
    
    calculate_Nq_();
    
    F_ = vector_dot(np_, Rpq_);
}


//virtual
dcomplex QuadKer_PlanVH_VolKer3_KerTyp2::specific_value_(const size_t ) const noexcept {
    
    return Nq_ * F_ * up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////


QuadKer_PlanVH_VolKer3_KerTyp3::QuadKer_PlanVH_VolKer3_KerTyp3() {
    up_green_func_.reset(new VolumeKernel3_KernelType3_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer3_KerTyp3::~QuadKer_PlanVH_VolKer3_KerTyp3() {
}

// virtual
void QuadKer_PlanVH_VolKer3_KerTyp3::precompute_rp_rq_dependent_data_()  noexcept {
    
    calculate_Np_();
    
    F_ = vector_dot(np_, Rpq_);
}


//virtual
dcomplex QuadKer_PlanVH_VolKer3_KerTyp3::specific_value_(const size_t ) const noexcept {
    
    return Np_ * F_ * up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////

QuadKer_PlanVH_VolKer3_KerTyp4::QuadKer_PlanVH_VolKer3_KerTyp4() {
    up_green_func_.reset(new VolumeKernel3_KernelType4_GreenFunc());
}

//virtual
QuadKer_PlanVH_VolKer3_KerTyp4::~QuadKer_PlanVH_VolKer3_KerTyp4() {
}

// virtual
void QuadKer_PlanVH_VolKer3_KerTyp4::precompute_rp_rq_dependent_data_()  noexcept {
}

//virtual
dcomplex QuadKer_PlanVH_VolKer3_KerTyp4::specific_value_(const size_t ) const noexcept {
    return up_green_func_->value();
}

///////////////////////////////////////////////////////////////////////////////

			
			
}  // End of the namespace Directfn

// End of the file


