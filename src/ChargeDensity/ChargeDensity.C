// File: ChargeDensity.C  Interface for a charge density
module;
#include <type_traits>
#include <atomic>
#include <cstddef>
export module qchem.ChargeDensity;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO
export import qchem.Symmetry.Spin;
import qchem.ScalarFunction;
import qchem.ChargeDensity.Types;

export namespace qchem::ChargeDensity
{

//! Monotonic logical-clock (Lamport) stamp for charge-density freshness.  Each new or mutated density
//! draws the next serial; the Hamiltonian term cache compares serials to ask "is this a *different*
//! density than the one I cached?".  It is NOT wall-clock time on purpose: a counter can never collide
//! (two densities born the same microsecond would), never runs backward (NTP/DST jumps would), and is
//! fully deterministic run-to-run / under any allocator (the property the whole determinism fix rests
//! on).  Per-T so the double and dcmplx lineages stay independent.  TRANSIENT: a density's serial is
//! runtime-only identity (like a pointer) -- never serialize it; a deserialized density re-stamps from
//! this counter.  Cross-run / value identity (a disk cache key) is a content hash, a separate concern.
template <class T> size_t NextDensityVersion()
{
    static std::atomic<size_t> theClock{0};
    return ++theClock;   // first serial handed out is 1; 0 is the reserved "no density yet" sentinel
}

//! Empty (non-polymorphic) stand-in for a PERIODIC density, which has no AO (auxiliary-basis) projection.
//! A density template inherits Fitting::ProjectedDensity_AO only on the finite (double) path; the periodic
//! (dcmplx) path gets this empty base, so it is NOT a ProjectedDensity_AO (FittedCD's cross-cast correctly
//! fails) and its object layout is unchanged.  The AO mirror of NoFourierDensity / FourierDensityBase: AO
//! projection (finite) and FT projection (periodic) are now both cross-cast capabilities, not forced bases.
struct NoProjectedDensity {};

//! ProjectedDensity_AO for the finite path (T=double), the empty base for the periodic path (T=dcmplx).
template <class T> using ProjectedDensityBase =
    std::conditional_t<std::is_same_v<T,double>, Fitting::ProjectedDensity_AO, NoProjectedDensity>;
//
//  These little interfaces allow us to invert a dependency with Hamiltonian Terms.
//  Templated on the matrix element type T (double for atoms/molecules; dcmplx for the
//  plane-wave lattice lineage where k-points make the blocks Hermitian-complex).  For T=double,
//  hmat_t<double> IS rsmat_t and tobs_t<double> IS obs_t, so the aliases below leave all existing
//  real code source- and binary-unchanged.
//
template <class T> class tDM_CD;   // forward (the contract clients name it)

template <class T> class tStatic_CC //Contract client for static Ham terms.
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&) const=0;
};

template <class T> class tDynamic_CC //Contract client for dynamic (CD dependent) Ham terms.
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&,const tDM_CD<T>*) const=0;
};

// Naming convention (mirrors rsmat_t/chmat_t in Common/Types.C): r* = <double>, c* = <dcmplx>.
// The bare names are TRANSITIONAL aliases to the r* version so existing real code is untouched;
// the full bare->r* rename across the codebase is pinned as a post-integration cleanup pass.
using rStatic_CC  = tStatic_CC<double>;   using cStatic_CC  = tStatic_CC<dcmplx>;
using rDynamic_CC = tDynamic_CC<double>;  using cDynamic_CC = tDynamic_CC<dcmplx>;
using Static_CC   = rStatic_CC;
using Dynamic_CC  = rDynamic_CC;

//----------------------------------------------------------------------------------
//
//  Charge density has a simple mandate:
//    1) Provide numerical evluation of ro(r).
//    2) Calculate the Coulomb self energy = sum ni <i(1)|Ro(2)/r12|i(1)> = sum Dab <a(1)|Ro(2)/r12|b(1)>
//    3) Calculate Vcoul(0) = <Ro(r)/r>.
//    4) Calculate the overlap   integrals  < ro(1)| b(1) > for some basis set b.
//    5) Calculate the repulsion integrals  < ro(1)/r12 | b(2) > for some basis set b.
//    6) Calculate the orbital repulsion integrals  < i(1) | ro(2)/r12 | j(1) > for orbitals i,j.
//    7) calculate the self repulsion = 1/2 <ro(1)|1/r12|ro(2)>
//
//  This is the interface for a charge density representation based on the density matrix.
//
template <class T> class tDM_CD
: public virtual ScalarFunction<double>
{
public:
    virtual double DM_Contract(const tStatic_CC<T>*) const=0; //Amounts to Integral(ro*V*d3r);
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const=0; //Integral(ro*V(ro)*d3r);

    virtual void   ReScale      (double factor              )      =0;  //Ro *= factor
    virtual void   MixIn        (const tDM_CD<T>&,double      )      =0;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tDM_CD<T>&            ) const=0;  //Convergence check.

    virtual double GetTotalCharge  () const=0;  // <ro>

    //! Logical-clock serial (see NextDensityVersion): distinct/mutated densities have distinct serials.
    //! TRANSIENT runtime identity -- not part of the persisted value (never serialize it).
    virtual size_t Version() const=0;

    // HF/fit-specific (real, Gaussian-basis) -- the dcmplx plane-wave density NA-asserts these.
    virtual void AccumulateDirect  (hmat_t<T>& Jab, const ohfbs_t*) const=0;
    virtual void AccumulateExchange(hmat_t<T>& Kab, const ohfbs_t*) const=0;

};

using rDM_CD = tDM_CD<double>;  
using cDM_CD = tDM_CD<dcmplx>;
using DM_CD  = rDM_CD;          // transitional bare alias

//---------------------------------------------------------------------------------------
//
//  Store spin up and spin down as a ChargeDensity
//  Generic: Could be fitted or exact.
//
class Polarized_CD
    : public virtual DM_CD
    , public virtual ProjectedDensityBase<double>   // finite/molecular: an AO-projectable density
{
public:
    virtual       DM_CD* GetChargeDensity(const Spin&)      =0;
    virtual const DM_CD* GetChargeDensity(const Spin&) const=0;

    virtual double DM_Contract(const Static_CC*) const;
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const;

    virtual double GetTotalCharge() const;  // <ro>
    virtual double GetTotalSpin  () const;  // No UT coverage// <up>-<down>

    // The spin children are mutated together (MixIn/ReScale below touch both), so either child's serial
    // tracks the polarized density's freshness; forward to Up.
    virtual size_t Version() const {return GetChargeDensity(Spin::Up)->Version();}

    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::FIT_CD_ABS*) const;
    virtual void AccumulateDirect  (rsmat_t& Jab, const ohfbs_t*) const;
    virtual void AccumulateExchange(rsmat_t& Kab, const ohfbs_t*) const;

    virtual void   ReScale      (double factor              )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //Convergence check.

    virtual double operator()(const rvec3_t&) const; // No UT coverage
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage
};

class SpinDensity : public virtual ScalarFunction<double>
{
public:
    SpinDensity(DM_CD* up,DM_CD* down);
    ~SpinDensity();
    virtual double operator()(const rvec3_t&) const; // No UT coverage
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage
private:
    DM_CD* itsSpinUpCD;
    DM_CD* itsSpinDownCD;
};

} //namespace