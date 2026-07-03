// File: ChargeDensity.C  Interface for a charge density
module;
#include <type_traits>
#include <cstddef>
#include <atomic>
#include <vector>
#include <map>
#include <string>
#include <cassert>
export module qchem.ChargeDensity;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO
export import qchem.Symmetry.Spin;
import qchem.ScalarFunction;
import qchem.ChargeDensity.Types;

export namespace qchem::ChargeDensity
{

//! Hand out the next TRANSIENT density-freshness serial -- the monotonic logical clock the dynamic-term
//! caches key on (Version()).  EVERY concrete density (IrrepCD, NumericCD, FourierSeedCD, ...) MUST
//! draw from this ONE counter: a serial that collides across density KINDS makes a dynamic term reuse a
//! stale cached matrix (the iter-0 seed Fock for the iter-1 working density), silently breaking the SCF --
//! see [[project_hamiltonian_dynamic_cache_bug]].  A SINGLE program-wide std::atomic (not per-T: the cache
//! only needs globally-distinct serials, and one counter makes EVERY density's Version() unique, real or
//! complex).  Serial 0 is the reserved "no density yet" sentinel; the first handed out is 1.
inline size_t NextDensityVersion()
{
    static std::atomic<size_t> theClock{0};
    return ++theClock;
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
template <class T> class tChargeDensity;   // forward (the DFT face the framework consumes)
template <class T> class tDM_CD;           // forward (the density-matrix face)

template <class T> class tStatic_CC //Contract client for static Ham terms.
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&) const=0;
};

template <class T> class tDynamic_CC //Contract client for dynamic (CD dependent) Ham terms.
{
public:
    // The density that defines V(rho) for the Fock build -- the DFT-only face (a fit has no matrix).
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*) const=0;
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
//  Base charge-density interface: everything that does NOT need a density matrix -- evaluate ro(r)
//  (ScalarFunction), the total charge, a transient freshness serial, and uniform scaling.  A
//  fitted/analytic density (e.g. the SAD seed, NumericCD) IS-A tChargeDensity but NOT a tDM_CD.
//  The DFT Fock build consumes densities through THIS face; the HF terms cross-cast to tDM_CD.
//
template <class T> class tChargeDensity
: public virtual ScalarFunction<double>
{
public:
    virtual double GetTotalCharge() const=0;       // <ro>
    virtual void   ReScale(double factor)     =0;  // Ro *= factor

    //! Monotonic logical-clock serial: distinct (or mutated) densities have distinct serials, so a cache
    //! can ask "is this a *different* density than the one I hold?".  TRANSIENT runtime identity (like a
    //! pointer) -- not part of the persisted value, never serialize it.  (Concrete densities stamp this
    //! from a per-T counter in IrrepCD's impl; composites/polarized forward to a child.)
    virtual size_t Version() const=0;
};

//  A charge density represented BY A DENSITY MATRIX: adds the matrix-only capabilities -- operator
//  contraction (DM_Contract), SCF density mixing (MixIn/GetChangeFrom), and the Hartree-Fock J/K
//  accumulators.  Densities with no matrix (fits) are tChargeDensity instead.
//
template <class T> class tDM_CD
: public virtual tChargeDensity<T>
{
public:
    virtual double DM_Contract(const tStatic_CC<T>*) const=0; //Amounts to Integral(ro*V*d3r);
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const=0; //Integral(ro*V(ro)*d3r);
    //! Energy contraction of THIS density against a per-irrep Fock-block map (keyed by BasisSetID):
    //! \f$\sum_i \mathrm{sum}(D_i \odot B_i^\mathsf{T}) = \mathrm{Tr}(D\,B)\f$.  The whole-system DUAL of
    //! AccumulateDirectAll -- lets a whole-system term (HF Coulomb/exchange) take its energy from its own
    //! cached blocks WITHOUT a per-irrep GetMatrix round-trip through DM_Contract.  Only composite/
    //! polarized/leaf densities implement it; the periodic path asserts out.
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const
    { assert(false && "DM_ContractBlocks: only a finite (Gaussian) density matrix implements this"); return 0.0; }

    virtual void   MixIn        (const tDM_CD<T>&,double      )      =0;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tDM_CD<T>&            ) const=0;  //Convergence check.

    // HF/fit-specific (real, Gaussian-basis) -- the dcmplx plane-wave density NA-asserts these.
    virtual void AccumulateDirect  (hmat_t<T>& Jab, const ohfbs_t*) const=0;
    virtual void AccumulateExchange(hmat_t<T>& Kab, const ohfbs_t*) const=0;

    //! Whole-system HF Coulomb exploiting the ERI4 bra-ket symmetry \f$J(i,j)=J(j,i)^\mathsf{T}\f$
    //! (doc/ERI4Rework.md \S4/\S5.4): one pass over each canonical ab-irrep pair feeds BOTH Fock blocks,
    //! so J(j,i) is never built or stored.  \a Jall is parallel to \a abBases (the Hamiltonian context's
    //! irrep bases) and pre-zeroed by the caller (the Vee term).  Only a composite / polarized density
    //! -- which spans every irrep block -- implements this; a lone leaf and the periodic path assert out.
    virtual void AccumulateDirectAll(std::vector<hmat_t<T>>& Jall, const std::vector<const ohfbs_t*>& abBases) const
    { assert(false && "AccumulateDirectAll: only a composite/polarized density spans all irrep blocks"); }
    //! Fused scatter of ONE canonical pair (this=irrep i, \a other=irrep j) into their two Fock blocks --
    //! the per-pair helper the composite's AccumulateDirectAll loop calls.  Only a leaf (IrrepCD) is a
    //! bra-ket pair partner; composite/polarized inherit this asserting default.
    virtual void AccumulateDirectBoth(hmat_t<T>& Ji, hmat_t<T>& Jj, const tDM_CD<T>& other) const
    { assert(false && "AccumulateDirectBoth: only a leaf (irrep) density is a bra-ket pair partner"); }
    //! Exchange analogues of the two above.  Exchange is SAME-SPIN, so the whole-system build is per single-
    //! spin density (Polarized_CD sums the two spin channels' K into the same blocks for the RHF term).
    virtual void AccumulateExchangeAll(std::vector<hmat_t<T>>& Kall, const std::vector<const ohfbs_t*>& abBases) const
    { assert(false && "AccumulateExchangeAll: only a composite/polarized density spans all irrep blocks"); }
    virtual void AccumulateExchangeBoth(hmat_t<T>& Ki, hmat_t<T>& Kj, const tDM_CD<T>& other) const
    { assert(false && "AccumulateExchangeBoth: only a leaf (irrep) density is a bra-ket pair partner"); }
};

using rChargeDensity = tChargeDensity<double>;  using cChargeDensity = tChargeDensity<dcmplx>;
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
    virtual double DM_ContractBlocks(const std::map<std::string,rsmat_t>&) const;   // sum both spins

    virtual double GetTotalCharge() const;  // <ro>
    virtual double GetTotalSpin  () const;  // No UT coverage// <up>-<down>

    // The spin children are mutated together (MixIn/ReScale below touch both), so either child's serial
    // tracks the polarized density's freshness; forward to Up.
    virtual size_t Version() const {return GetChargeDensity(Spin::Up)->Version();}

    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::FIT_CD_ABS*) const;
    virtual void AccumulateDirect  (rsmat_t& Jab, const ohfbs_t*) const;
    virtual void AccumulateExchange(rsmat_t& Kab, const ohfbs_t*) const;
    virtual void AccumulateDirectAll  (std::vector<rsmat_t>& Jall, const std::vector<const ohfbs_t*>& abBases) const;  // sum both spins (Coulomb)
    virtual void AccumulateExchangeAll(std::vector<rsmat_t>& Kall, const std::vector<const ohfbs_t*>& abBases) const;  // sum both spins (RHF exchange)

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