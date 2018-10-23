#include <iostream>
#include <iomanip>
#include "directfn_contour.h"

using  std::cout;
using  std::cerr;
using  std::endl;
using  std::setw;
using  std::setprecision;

namespace Directfn {

///////////////////////////////////////////////////////////////////////////////

SingularContour3xn::SingularContour3xn():
// The data_ array is too long and here is not convenient
// to init all the 21 zeros. It is done bellow in the loop.
r1_plus_r3_{0.0, 0.0, 0.0},
r4_plus_r6_{0.0, 0.0, 0.0},
r3_plus_r6_{0.0, 0.0, 0.0} {

    // Zero initialization
    for (size_t i = 0; i < N3_() * 7; ++i) {
        data_[i] = 0.0;
    }
    length_ = 0;
}

SingularContour3xn::~SingularContour3xn() {
}

/*! Setup contour for Self-Term interaction. */
void SingularContour3xn::set_points(const double * p1, const double * p2, const double * p3) {

    for (size_t i = 0; i < N3_(); ++i) {
        data_[0 * N3_() + i] = p1[i];
        data_[1 * N3_() + i] = p2[i];
        data_[2 * N3_() + i] = p3[i];
    }
    length_ = 3;
}

void SingularContour3xn::set_points(const double * p1, const double * p2, const double * p3, const double * p4) {

    for (size_t i = 0; i < N3_(); ++i) {
        data_[0 * N3_() + i] = p1[i];
        data_[1 * N3_() + i] = p2[i];
        data_[2 * N3_() + i] = p3[i];
        data_[3 * N3_() + i] = p4[i];
    }
    length_ = 4;

    update_faces_centers_st_();
}

void SingularContour3xn::set_points(const double * p1, const double * p2, const double * p3,
                                    const double * p4, const double * p5) {

    for (size_t i = 0; i < N3_(); ++i) {
        data_[0 * N3_() + i] = p1[i];
        data_[1 * N3_() + i] = p2[i];
        data_[2 * N3_() + i] = p3[i];
        data_[3 * N3_() + i] = p4[i];
        data_[4 * N3_() + i] = p5[i];
    }
    length_ = 5;
}

void SingularContour3xn::set_points(const double * p1, const double * p2, const double * p3, const double * p4,
                                    const double * p5, const double * p6) {
    this->set_points(p1, p2, p3, p4);

    for (size_t i = 0; i < N3_(); ++i) {
        data_[4 * N3_() + i] = p5[i];
        data_[5 * N3_() + i] = p6[i];
    }
    length_ = 6;

    update_faces_centers_ea_();
}

void SingularContour3xn::set_points(const double * p1, const double * p2, const double * p3, const double * p4,
                                    const double * p5, const double * p6, const double * p7) {

    this->set_points(p1, p2, p3, p4, p5, p6);

    for (size_t i = 0; i < N3_(); ++i) {
        data_[6 * N3_() + i] = p7[i];
    }
    length_ = 7;

    update_faces_centers_va_();
}

void SingularContour3xn::set_points(const double * p1, const double * p2, const double * p3,
                                    const double * p4, const double * p5, const double * p6,
                                    const double * p7, const double * p8, const double * p9) {

    for (size_t i = 0; i < N3_(); ++i) {

        data_[0 * N3_() + i] = p1[i];
        data_[1 * N3_() + i] = p2[i];
        data_[2 * N3_() + i] = p3[i];
        data_[3 * N3_() + i] = p4[i];
        data_[4 * N3_() + i] = p5[i];
        data_[5 * N3_() + i] = p6[i];
        data_[6 * N3_() + i] = p7[i];
        data_[7 * N3_() + i] = p8[i];
        data_[8 * N3_() + i] = p9[i];
    }
    length_ = 9;
}

void SingularContour3xn::set_points(const double * p1,  const double * p2,  const double * p3,
                                    const double * p4,  const double * p5,  const double * p6,
                                    const double * p7,  const double * p8,  const double * p9,
                                    const double * p10, const double * p11, const double * p12,
                                    const double * p13, const double * p14, const double * p15) {

    this->set_points(p1, p2, p3, p4, p5, p6, p7, p8, p9);

    for (size_t i = 0; i < N3_(); ++i) {

        data_[ 9 * N3_() + i] = p10[i];
        data_[10 * N3_() + i] = p11[i];
        data_[11 * N3_() + i] = p12[i];
        data_[12 * N3_() + i] = p13[i];
        data_[13 * N3_() + i] = p14[i];
        data_[14 * N3_() + i] = p15[i];
    }
    length_ = 15;
}

void  SingularContour3xn::set_points(const double * p1,  const double * p2,  const double * p3,
                                     const double * p4,  const double * p5,  const double * p6,
                                     const double * p7,  const double * p8,  const double * p9,
                                     const double * p10, const double * p11, const double * p12,
                                     const double * p13, const double * p14, const double * p15,
                                     const double * p16, const double * p17) {

    this->set_points(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15);

    for (size_t i = 0; i < N3_(); ++i) {

        data_[15 * N3_() + i] = p16[i];
        data_[16 * N3_() + i] = p17[i];
    }

    length_ = 17;
}

size_t  SingularContour3xn::length() const noexcept {
    return length_;
}

const double * SingularContour3xn::operator() (const size_t index) const noexcept {

    if (index >= 18) {
        cerr << "Error 1154-14122014" << endl;
        // TO DO:
        // change to message logger...
    }
    return &data_[N3_() * index];
}

double SingularContour3xn::x(const size_t index) const noexcept {
    return data_[N3_() * index + 0];
}

double SingularContour3xn::y(const size_t index) const noexcept {
    return data_[N3_() * index + 1];
}

double SingularContour3xn::z(const size_t index) const noexcept {
    return data_[N3_() * index + 2];
}

/*! Return pointer to the r1_plus_r3_ [3] - array. */
const double * SingularContour3xn::r1_plus_r3_half() const noexcept {
    return r1_plus_r3_;
}

const double * SingularContour3xn::r4_plus_r6_half() const noexcept {
    return r4_plus_r6_;
}

const double * SingularContour3xn::r3_plus_r6_half() const noexcept {
    return r3_plus_r6_;
}

void SingularContour3xn::debug_print() const noexcept {

    for (size_t k = 0; k < length(); ++k) { cout << setw(10) << x(k);} cout << endl;
    for (size_t k = 0; k < length(); ++k) { cout << setw(10) << y(k);} cout << endl;
    for (size_t k = 0; k < length(); ++k) { cout << setw(10) << z(k);} cout << endl;
}

void SingularContour3xn::update_faces_centers_st_() noexcept {

    const double * const r1 = &data_[N3_() * 0 + 0];
    const double * const r3 = &data_[N3_() * 2 + 0];

    for (size_t i = 0; i < N3_(); ++i) {
        r1_plus_r3_[i] = 0.5 * (r1[i] + r3[i]);
    }
}

void SingularContour3xn::update_faces_centers_ea_() noexcept {

    const double * const r1 = &data_[N3_() * 0 + 0];
    const double * const r3 = &data_[N3_() * 2 + 0];
    const double * const r4 = &data_[N3_() * 3 + 0];
    const double * const r6 = &data_[N3_() * 5 + 0];

    for (size_t i = 0; i < N3_(); ++i) {
        r1_plus_r3_[i] = 0.5 * (r1[i] + r3[i]);
        r4_plus_r6_[i] = 0.5 * (r4[i] + r6[i]);
    }
}

void SingularContour3xn::update_faces_centers_va_() noexcept {

    const double * const r1 = &data_[N3_() * 0 + 0];
    const double * const r3 = &data_[N3_() * 2 + 0];
    const double * const r6 = &data_[N3_() * 5 + 0];

    for (size_t i = 0; i < N3_(); ++i) {
        r1_plus_r3_[i] = 0.5 * (r1[i] + r3[i]);
        r3_plus_r6_[i] = 0.5 * (r3[i] + r6[i]);
    }
}

///////////////////////////////////////////////////////////////////////////////

} // End of the namespace Directfn

// End of the file

