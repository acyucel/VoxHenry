#include "directfn_triag.h"

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_tri_st_plan(const double r1[3], const double r2[3], const double r3[3],
                         const size_t N1, const size_t N2,
                         const size_t N3, const size_t N4,
                         const double k0wn, dcomplex * const cp_data) {

    // In cp_data the size is supposed to be setup properly

    // Setup the integration contour verices
    static SingularContour3xn cntr3p;
    cntr3p.set_points(r1, r2, r3);

    static std::unique_ptr<Triangular_ST<ParticularKernel>> up_tri_st(new Triangular_ST<ParticularKernel>());
    up_tri_st->set_wavenumber(k0wn);
    up_tri_st->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_tri_st->set(cntr3p);

    // Computational routine
    up_tri_st->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_tri_st->kernel_size()]);
    up_tri_st->copy_Iss_array_values_to(cp_data);

    // Successfull return
    return 0;
}

// Instantiation
template int directfn_tri_st_plan<TriangularKernel_Constant_ST>(const double[3] , const double[3] , const double[3] ,
                                                                const size_t , const size_t ,
                                                                const size_t , const size_t ,
                                                                const double , dcomplex * const );

template int directfn_tri_st_plan<TriangularKernel_RWG_WS>(const double[3] , const double[3] , const double[3] ,
                                                           const size_t , const size_t ,
                                                           const size_t , const size_t ,
                                                           const double , dcomplex * const );

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_tri_ea_plan(const double r1[3], const double r2[3],
                         const double r3[3], const double r4[3],
                         const size_t N1, const size_t N2,
                         const size_t N3, const size_t N4,
                         const double k0wn, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr4p;
    cntr4p.set_points(r1, r2, r3, r4);

    static std::unique_ptr<Triangular_EA<ParticularKernel>> up_tri_ea(new Triangular_EA<ParticularKernel>());
    up_tri_ea->set_wavenumber(k0wn);
    up_tri_ea->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_tri_ea->set(cntr4p);

    // Computational routine
    up_tri_ea->calc_Iss();

    // Copy computed values to the buffer
    //up_data.reset(new dcomplex[up_tri_ea->kernel_size()]);
    up_tri_ea->copy_Iss_array_values_to(cp_data);

    // Successfull return
    return 0;
}


// Instantiation
template int directfn_tri_ea_plan<TriangularKernel_Constant_EA>(const double[3], const double[3],
                                                                const double[3], const double[3],
                                                                const size_t, const size_t,
                                                                const size_t, const size_t,
                                                                const double, dcomplex * const cp_data);

template int directfn_tri_ea_plan<TriangularKernel_RWG_WS>(const double[3], const double[3],
                                                           const double[3], const double[3],
                                                           const size_t, const size_t,
                                                           const size_t, const size_t,
                                                           const double, dcomplex * const cp_data);

template int directfn_tri_ea_plan<TriangularKernel_RWG_SS>(const double[3], const double[3],
                                                           const double[3], const double[3],
                                                           const size_t, const size_t,
                                                           const size_t, const size_t,
                                                           const double, dcomplex * const cp_data);

template int directfn_tri_ea_plan<TriangularKernel_nxRWG_SS>(const double[3], const double[3],
                                                             const double[3], const double[3],
                                                             const size_t, const size_t,
                                                             const size_t, const size_t,
                                                             const double, dcomplex * const cp_data);

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
int directfn_tri_va_plan(const double r1[3], const double r2[3],
                         const double r3[3], const double r4[3],
                         const double r5[3],
                         const size_t N1, const size_t N2,
                         const size_t N3, const size_t N4,
                         const double k0wn, dcomplex * const cp_data) {

    // Setup the integration contour verices
    static SingularContour3xn cntr5p;
    cntr5p.set_points(r1, r2, r3, r4, r5);

    static std::unique_ptr<Triangular_VA<ParticularKernel>> up_tri_va(new Triangular_VA<ParticularKernel>());
    up_tri_va->set_wavenumber(k0wn);
    up_tri_va->set_Gaussian_orders_4(N1, N2, N3, N4);
    up_tri_va->set(cntr5p);

    // Computational routine
    up_tri_va->calc_Iss();

    // Copy computed values to the buffer
    // up_data.reset(new dcomplex[up_tri_va->kernel_size()]);
    up_tri_va->copy_Iss_array_values_to(cp_data);

    // Successfull return the size of the array;
    return 0; //up_tri_va->kernel_size();
}


// Instantiation
template int directfn_tri_va_plan<TriangularKernel_Constant_VA>(const double[3], const double[3],
                                                                const double[3], const double[3],
                                                                const double[3],
                                                                const size_t, const size_t,
                                                                const size_t, const size_t,
                                                                const double, dcomplex * const );

template int directfn_tri_va_plan<TriangularKernel_RWG_WS>(const double[3], const double[3],
                                                           const double[3], const double[3],
                                                           const double[3],
                                                           const size_t, const size_t,
                                                           const size_t, const size_t,
                                                           const double, dcomplex * const );

template int directfn_tri_va_plan<TriangularKernel_RWG_SS>(const double[3], const double[3],
                                                           const double[3], const double[3],
                                                           const double[3],
                                                           const size_t, const size_t,
                                                           const size_t, const size_t,
                                                           const double, dcomplex * const );

template int directfn_tri_va_plan<TriangularKernel_nxRWG_SS>(const double[3], const double[3],
                                                             const double[3], const double[3],
                                                             const double[3],
                                                             const size_t, const size_t,
                                                             const size_t, const size_t,
                                                             const double, dcomplex * const );

///////////////////////////////////////////////////////////////////////////////
}  // End of the namespace Directfn

// End of the file

