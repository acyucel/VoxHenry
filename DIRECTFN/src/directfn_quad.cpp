#include "directfn_quad.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_st_plan(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0wn, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr4p;
    cntr4p.set_points(r1, r2, r3, r4);

    static std::unique_ptr<Quadrilateral_ST<ParticularKernel>> up_quad_st(new Quadrilateral_ST<ParticularKernel>());
    up_quad_st->set_wavenumber(k0wn);
    up_quad_st->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_st->set(cntr4p);

    // Computational routine
    up_quad_st->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_quad_st->kernel_size()]);
    up_quad_st->copy_Iss_array_values_to(cp_data);

    // Successfull return
    return 0;
}


// Instantiation
template int directfn_quad_st_plan<QuadrilateralKernel_PlanarScalar>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , dcomplex * const );

template int directfn_quad_st_plan<QuadrilateralKernel_PlanarVectorWS>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , dcomplex * const );

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_ea_plan(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0wn, dcomplex * const cp_data){

    // Setup the integration contour verices
    static SingularContour3xn cntr6p;
    cntr6p.set_points(r1, r2, r3, r4, r5, r6);

    static std::unique_ptr<Quadrilateral_EA<ParticularKernel>> up_quad_ea(new Quadrilateral_EA<ParticularKernel>());
    up_quad_ea->set_wavenumber(k0wn);
    up_quad_ea->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_ea->set(cntr6p);

    // Computational routine
    up_quad_ea->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_quad_ea->kernel_size()]);
    up_quad_ea->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0;  // up_quad_ea->kernel_size();
}


// Instantiation
template int directfn_quad_ea_plan<QuadrilateralKernel_PlanarScalar>(const double [3], const double [3],
                                                                     const double [3], const double [3],
                                                                     const double [3], const double [3],
                                                                     const size_t , const size_t ,
                                                                     const size_t , const size_t ,
                                                                     const double , dcomplex * const );

template int directfn_quad_ea_plan<QuadrilateralKernel_PlanarVectorWS>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , dcomplex * const );

template int directfn_quad_ea_plan<QuadrilateralKernel_PlanarVectorSS>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , dcomplex * const );

/////////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_va_plan(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const double r7[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0wn, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr7p;
    cntr7p.set_points(r1, r2, r3, r4, r5, r6, r7);

    static std::unique_ptr<Quadrilateral_VA<ParticularKernel>> up_quad_va(new Quadrilateral_VA<ParticularKernel>());
    up_quad_va->set_wavenumber(k0wn);
    up_quad_va->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_va->set(cntr7p);

    // Computational routine
    up_quad_va->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_quad_va->kernel_size()]);
    up_quad_va->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0;  // up_quad_va->kernel_size();
}


// Instantiation
template int directfn_quad_va_plan<QuadrilateralKernel_PlanarScalar>(const double [3], const double [3],
                                                                     const double [3], const double [3],
                                                                     const double [3], const double [3],
                                                                     const double [3],
                                                                     const size_t , const size_t ,
                                                                     const size_t , const size_t ,
                                                                     const double , dcomplex * const );

template int directfn_quad_va_plan<QuadrilateralKernel_PlanarVectorWS>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , dcomplex * const );

template int directfn_quad_va_plan<QuadrilateralKernel_PlanarVectorSS>(const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3], const double [3],
                                                                       const double [3],
                                                                       const size_t , const size_t ,
                                                                       const size_t , const size_t ,
                                                                       const double , dcomplex * const );

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_st_curv(const double r1[3], const double r2[3],
                          const double r3[3], const double r4[3],
                          const double r5[3], const double r6[3],
                          const double r7[3], const double r8[3],
                          const double r9[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0wn, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr9p;
    cntr9p.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9);

    static std::unique_ptr<Quadrilateral_ST<ParticularKernel>> up_quad_st(new Quadrilateral_ST<ParticularKernel>());
    up_quad_st->set_wavenumber(k0wn);
    up_quad_st->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_st->set(cntr9p);

    // Computational routine
    up_quad_st->calc_Iss();

    // Copy computed values to the buffer
//    up_data.reset(new dcomplex[up_quad_st->kernel_size()]);
    up_quad_st->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0;  //up_quad_st->kernel_size();
}


// Instantiation
template int directfn_quad_st_curv<QuadrilateralKernel_CurvilinearScalar>(const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3],
                                                                          const size_t , const size_t ,
                                                                          const size_t , const size_t ,
                                                                          const double , dcomplex * const);

template int directfn_quad_st_curv<QuadrilateralKernel_CurvilinearVectorWS>(const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3],
                                                                            const size_t , const size_t ,
                                                                            const size_t , const size_t ,
                                                                            const double , dcomplex * const);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_ea_curv(const double r1[3],  const double r2[3],
                          const double r3[3],  const double r4[3],
                          const double r5[3],  const double r6[3],
                          const double r7[3],  const double r8[3],
                          const double r9[3],  const double r10[3],
                          const double r11[3], const double r12[3],
                          const double r13[3], const double r14[3],
                          const double r15[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr15p;
    cntr15p.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15);

    static std::unique_ptr<Quadrilateral_EA<ParticularKernel>> up_quad_ea(new Quadrilateral_EA<ParticularKernel>());
    up_quad_ea->set_wavenumber(k0);
    up_quad_ea->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_ea->set(cntr15p);

    // Computational routine
    up_quad_ea->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_quad_ea->kernel_size()]);
    up_quad_ea->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0; //up_quad_ea->kernel_size();
}

// Instantiation
template int directfn_quad_ea_curv<QuadrilateralKernel_CurvilinearScalar>(const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3],
                                                                          const size_t , const size_t ,
                                                                          const size_t , const size_t ,
                                                                          const double , dcomplex * const );

template int directfn_quad_ea_curv<QuadrilateralKernel_CurvilinearVectorWS>(const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3],
                                                                            const size_t , const size_t ,
                                                                            const size_t , const size_t ,
                                                                            const double , dcomplex * const );

template int directfn_quad_ea_curv<QuadrilateralKernel_CurvilinearVectorSS>(const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3],
                                                                            const size_t , const size_t ,
                                                                            const size_t , const size_t ,
                                                                            const double , dcomplex * const );

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_quad_va_curv(const double r1[3],  const double r2[3],
                          const double r3[3],  const double r4[3],
                          const double r5[3],  const double r6[3],
                          const double r7[3],  const double r8[3],
                          const double r9[3],  const double r10[3],
                          const double r11[3], const double r12[3],
                          const double r13[3], const double r14[3],
                          const double r15[3], const double r16[3],
                          const double r17[3],
                          const size_t N1, const size_t N2,
                          const size_t N3, const size_t N4,
                          const double k0, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr17p;
    cntr17p.set_points(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17);

    static std::unique_ptr<Quadrilateral_VA<ParticularKernel>> up_quad_va(new Quadrilateral_VA<ParticularKernel>());
    up_quad_va->set_wavenumber(k0);
    up_quad_va->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_quad_va->set(cntr17p);

    // Computational routine
    up_quad_va->calc_Iss();

    // Copy computed values to the buffer
//    up_data.reset(new dcomplex[up_quad_va->kernel_size()]);
    up_quad_va->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0; // up_quad_va->kernel_size();
}


// Instantiation
template int directfn_quad_va_curv<QuadrilateralKernel_CurvilinearScalar>(const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3], const double [3],
                                                                          const double [3],
                                                                          const size_t , const size_t ,
                                                                          const size_t , const size_t ,
                                                                          const double , dcomplex * const);

template int directfn_quad_va_curv<QuadrilateralKernel_CurvilinearVectorWS>(const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3],
                                                                            const size_t , const size_t ,
                                                                            const size_t , const size_t ,
                                                                            const double , dcomplex * const);

template int directfn_quad_va_curv<QuadrilateralKernel_CurvilinearVectorSS>(const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3], const double [3],
                                                                            const double [3],
                                                                            const size_t , const size_t ,
                                                                            const size_t , const size_t ,
                                                                            const double , dcomplex * const);

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

// End of the file

