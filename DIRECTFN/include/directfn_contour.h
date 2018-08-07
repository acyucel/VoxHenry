#ifndef _DIRECTFN_SINGULAR_CONTOUR_H_
#define _DIRECTFN_SINGULAR_CONTOUR_H_

#include <cstddef>

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

class SingularContour3xn final {
public:
    SingularContour3xn();
    ~SingularContour3xn();

    SingularContour3xn(const SingularContour3xn & ) = default;
    SingularContour3xn & operator = (const SingularContour3xn & ) = default;

    SingularContour3xn(SingularContour3xn && ) = delete;
    SingularContour3xn & operator = (SingularContour3xn && ) = delete;

    /*! Setup contour for Self-Term interaction. */
    void  set_points(const double * p1, const double * p2, const double * p3);

    /*! Setup contour for Self-Term interaction. */
    void  set_points(const double * p1, const double * p2, const double * p3, const double * p4);

    void  set_points(const double * p1, const double * p2, const double * p3,
                     const double * p4, const double * p5);

    /*! Setup contour for Edge-adjacent interaction. */
    void  set_points(const double * p1, const double * p2, const double * p3, const double * p4,
                     const double * p5, const double * p6);

    /*! Setup contour for Vertex adjacent interaction. */
    void  set_points(const double * p1, const double * p2, const double * p3, const double * p4,
                     const double * p5, const double * p6, const double * p7);

    /*! Setup contour for Vertex adjacent interaction for Curvilinear elements. */
    void  set_points(const double * p1, const double * p2, const double * p3,
                     const double * p4, const double * p5, const double * p6,
                     const double * p7, const double * p8, const double * p9);

    /*! Curvilinear EA */
    void  set_points(const double * p1,  const double * p2,  const double * p3,
                     const double * p4,  const double * p5,  const double * p6,
                     const double * p7,  const double * p8,  const double * p9,
                     const double * p10, const double * p11, const double * p12,
                     const double * p13, const double * p14, const double * p15);

    /*! Curvilinear VA */
    void  set_points(const double * p1,  const double * p2,  const double * p3,
                     const double * p4,  const double * p5,  const double * p6,
                     const double * p7,  const double * p8,  const double * p9,
                     const double * p10, const double * p11, const double * p12,
                     const double * p13, const double * p14, const double * p15,
                     const double * p16, const double * p17);

    size_t  length() const noexcept;

    const double * operator() (const size_t index) const noexcept;

    double x(const size_t index) const noexcept;
    double y(const size_t index) const noexcept;
    double z(const size_t index) const noexcept;

    /*! \return pointer to the r1_plus_r3 buffer.
     *  call these 3 carefully - can be used only for
     *  st or va or ea after set_points */
    const double * r1_plus_r3_half() const noexcept;
    const double * r4_plus_r6_half() const noexcept;
    const double * r3_plus_r6_half() const noexcept;

    void  debug_print() const noexcept;

private:
    static constexpr size_t N3_() noexcept {return 3;}

    /*! 17 = 9 + 8 for VA Curvilinear elements  */
    double  data_[3 * 17];
    size_t  length_;
    /*! This var is needed for ST quadrilateral integration over faces of the cube.
     *  For the other casesmayve the center of the face should be changed.  */
    double  r1_plus_r3_[3];
    double  r4_plus_r6_[3];
    double  r3_plus_r6_[3];

    void  update_faces_centers_st_() noexcept;
    void  update_faces_centers_ea_() noexcept;
    void  update_faces_centers_va_() noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif  // _DIRECTFN_SINGULAR_CONTOUR_H_

// End of the file

