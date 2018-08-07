#include "directfn_quad_voxhenry.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_st_plan_voxhenry(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0wn, const double nq[3], const double np[3],
                          const int lp, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr4p;
    cntr4p.set_points(r1, r2, r3, r4);

    static std::unique_ptr<Quadrilateral_ST_VH<ParticularKernel>> up_quad_st(new Quadrilateral_ST_VH<ParticularKernel>());
    up_quad_st->set_wavenumber(k0wn);
    up_quad_st->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_st->set(cntr4p);
    up_quad_st->set_lp(lp);
    up_quad_st->set_nq_np(nq, np);

    // Computational routine
    up_quad_st->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_quad_st->kernel_size()]);
    up_quad_st->copy_Iss_array_values_to(cp_data);

    // Successfull return
    return 0;
}


// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp1>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp2>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );


// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp3>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );


// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp4>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp1>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp2>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp3>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_st_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp4>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );


///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_ea_plan_voxhenry(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0wn, const double nq[3], const double np[3],
                          const int lp, dcomplex * const cp_data){

    // Setup the integration contour verices
    static SingularContour3xn cntr6p;
    cntr6p.set_points(r1, r2, r3, r4, r5, r6);

    static std::unique_ptr<Quadrilateral_EA_VH<ParticularKernel>> up_quad_ea(new Quadrilateral_EA_VH<ParticularKernel>());
    up_quad_ea->set_wavenumber(k0wn);
    up_quad_ea->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_ea->set(cntr6p);
    up_quad_ea->set_lp(lp);
    up_quad_ea->set_nq_np(nq, np);
        
    // Computational routine
    up_quad_ea->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_quad_ea->kernel_size()]);
    up_quad_ea->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0;  // up_quad_ea->kernel_size();
}


// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp1>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp2>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp3>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp4>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp1>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp2>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp3>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_ea_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp4>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );



/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_va_plan_voxhenry(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const double r7[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0wn, const double nq[3], const double np[3],
                          const int lp, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr7p;
    cntr7p.set_points(r1, r2, r3, r4, r5, r6, r7);

    static std::unique_ptr<Quadrilateral_VA_VH<ParticularKernel>> up_quad_va(new Quadrilateral_VA_VH<ParticularKernel>());
    up_quad_va->set_wavenumber(k0wn);
    up_quad_va->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_va->set(cntr7p);
    up_quad_va->set_lp(lp);
    up_quad_va->set_nq_np(nq, np);
    
    // Computational routine
    up_quad_va->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_quad_va->kernel_size()]);
    up_quad_va->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0;  // up_quad_va->kernel_size();
}


// Instantiation
// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp1>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp2>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp3>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer2_KerTyp4>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp1>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp2>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp3>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );

// Instantiation
template int directfn_quad_va_plan_voxhenry<QuadKer_PlanVH_VolKer3_KerTyp4>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , const double [3], const double [3],
                                                                       const int , dcomplex * const );



}  // End of the namespace Directfn

// End of the file

