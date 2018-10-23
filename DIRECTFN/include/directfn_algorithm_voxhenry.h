#ifndef _EDGEADJACENT_ANYLATERAL_VH_H_
#define _EDGEADJACENT_ANYLATERAL_VH_H_

#include "directfn_algorithm_ea.h"
#include "directfn_algorithm_st.h"
#include "directfn_algorithm_va.h"

namespace Directfn {

// should implement in "directfn_interface.h" to avoid declaring the same function three times,
// but want to avoid modifying the DIRECTFN sources as well, so just deriving classes from it

template <typename ParticularKernel>
class Quadrilateral_EA_VH : public Quadrilateral_EA<ParticularKernel> {
public:
    Quadrilateral_EA_VH();
    ~Quadrilateral_EA_VH();

    Quadrilateral_EA_VH(const Quadrilateral_EA_VH & ) = delete;
    Quadrilateral_EA_VH(Quadrilateral_EA_VH && ) = delete;
    Quadrilateral_EA_VH & operator = (const Quadrilateral_EA_VH & ) = delete;
    Quadrilateral_EA_VH & operator = (Quadrilateral_EA_VH && ) = delete;

    void set_lp(const int lp) noexcept;
    void set_nq_np(const double nq[3], const double np[3]) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class Quadrilateral_ST_VH : public Quadrilateral_ST<ParticularKernel> {
public:
    Quadrilateral_ST_VH();
    ~Quadrilateral_ST_VH();

    Quadrilateral_ST_VH(const Quadrilateral_ST_VH & ) = delete;
    Quadrilateral_ST_VH(Quadrilateral_ST_VH && ) = delete;
    Quadrilateral_ST_VH & operator = (const Quadrilateral_ST_VH & ) = delete;
    Quadrilateral_ST_VH & operator = (Quadrilateral_ST_VH && ) = delete;

    void set_lp(const int lp) noexcept;
    void set_nq_np(const double nq[3], const double np[3]) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

template <typename ParticularKernel>
class Quadrilateral_VA_VH : public Quadrilateral_VA<ParticularKernel> {
public:
    Quadrilateral_VA_VH();
    ~Quadrilateral_VA_VH();

    Quadrilateral_VA_VH(const Quadrilateral_VA_VH & ) = delete;
    Quadrilateral_VA_VH(Quadrilateral_VA_VH && ) = delete;
    Quadrilateral_VA_VH & operator = (const Quadrilateral_VA_VH & ) = delete;
    Quadrilateral_VA_VH & operator = (Quadrilateral_VA_VH && ) = delete;

    void set_lp(const int lp) noexcept;
    void set_nq_np(const double nq[3], const double np[3]) noexcept;
};

///////////////////////////////////////////////////////////////////////////////

}  // End of the namespace Directfn

#endif   // _EDGEADJACENT_ANYLATERAL_VH_H_

// End of the file

