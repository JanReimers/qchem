// File: ChargeDensity.C  Interface for a charge density
module;
#include <type_traits>
#include <cstddef>
#include <atomic>
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cassert>
export module qchem.ChargeDensity;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO
export import qchem.Symmetry.Spin;
export import qchem.ChargeDensity.FourierDensity;   // FourierDensityBase<T> (tPolarized_CD's periodic face)
import qchem.ScalarFunction;
import qchem.ChargeDensity.Types;

export namespace qchem::ChargeDensity
{

//! Hand out the next TRANSIENT density-freshness serial -- the monotonic logical clock the dynamic-term
//! caches key on (Version()).  EVERY concrete density (IrrepCD, NumericCD, SeedCD, ...) MUST
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

//! A "lineage" is one charge density as it EVOLVES across an SCF run: iteration k's density is superseded by
//! iteration k+1's.  Consecutive densities are UNRELATED objects (each iteration builds a fresh one from the
//! new orbitals), so "which one is current?" cannot be answered by a density alone -- we track it explicitly.
//! This shared object holds one number: the Version() of whichever density is the live head.  A density is
//! active iff its Version() still equals that head (see tLineageTracked); an older, superseded density (a
//! kept itsOldCD, or a stale copy) has Version() < head and reports isActive()==false, so handing it to the
//! Hamiltonian trips an assert instead of silently building a Fock from the wrong density.  Head 0 = none yet.
class Lineage
{
    std::atomic<size_t> itsHead{0};
public:
    void   MakeHead(size_t v) {itsHead.store(v);}   //!< \a v is a NextDensityVersion serial (globally monotonic)
    size_t Head() const       {return itsHead.load();}
};
using LineagePtr = std::shared_ptr<Lineage>;

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
//  hmat_t<double> IS rsmat_t and tobs_t<double> IS robs_t, so the aliases below leave all existing
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
using rStatic_CC  = tStatic_CC<double>;   using cStatic_CC  = tStatic_CC<dcmplx>;
using rDynamic_CC = tDynamic_CC<double>;  using cDynamic_CC = tDynamic_CC<dcmplx>;

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

    //! \brief \f$\rho\f$ at MANY points -- the batched ScalarFunction.  Default: the pointwise loop
    //! (correct for any density); a concrete density overrides when it can vectorize the batch
    //! (FourierMixCD's factorized phase tables; the DM densities have the richer, table-driven
    //! DM_RhoAtPoints).  The mesh-XC quadrature engine's entry for non-DM densities.
    virtual rvec_t EvalBatch(const rvec3vec_t& r) const
    {
        rvec_t ro(r.size());
        for (size_t g=0; g<r.size(); g++) ro[g]=(*this)(r[g]);
        return ro;
    }

    //! Monotonic logical-clock serial: distinct (or mutated) densities have distinct serials, so a cache
    //! can ask "is this a *different* density than the one I hold?".  TRANSIENT runtime identity (like a
    //! pointer) -- not part of the persisted value, never serialize it.  (Concrete densities stamp this
    //! from a per-T counter in IrrepCD's impl; composites/polarized forward to a child.)
    virtual size_t Version() const=0;

    //! Layer-2 SCF-lineage check: am I the live head of my lineage, or a superseded density?  (See Lineage.)
    //! The default is "yes" -- a density that never joined a lineage (a seed, a one-off fit) is never a stale
    //! working density, so it is trivially active.  The top-level SCF working densities override this via
    //! tLineageTracked.  Consumers (the Hamiltonian) assert this before building a Fock from the density.
    virtual bool isActive() const {return true;}
    //! Become the current head of \a l (called by the SCFIterator on each new working density).  No-op for
    //! untracked densities.
    virtual void JoinLineage(const LineagePtr&) {}
};

//! Mixin that makes a versioned density track its SCF lineage head (see Lineage).  Mixed into the top-level
//! densities the SCFIterator hands the Hamiltonian (tComposite_CD, Polarized_CD).  A tracked density starts
//! with NO lineage (isActive() still true); JoinLineage makes it the head; MixIn/ReScale call AdvanceHead so
//! the mutated density stays the head (its Version() moved, so the head must move with it).
template <class T> class tLineageTracked : public virtual tChargeDensity<T>
{
    LineagePtr itsLineage;   //!< null until this density joins an SCF lineage
protected:
    void AdvanceHead() {if(itsLineage) itsLineage->MakeHead(this->Version());}   //!< after an in-place mutation
public:
    virtual bool isActive() const override {return !itsLineage || this->Version()==itsLineage->Head();}
    virtual void JoinLineage(const LineagePtr& l) override {itsLineage=l; l->MakeHead(this->Version());}
};

//! \brief A charge density that can be MIXED IN PLACE -- the whole of what an SCF density mixer needs of
//! its subject, and nothing else.
//!
//! ISP (user, 2026-08-08).  \c LinearMixer uses exactly \c GetChangeFrom, \c MixIn and \c GetTotalCharge,
//! yet it used to take a \c tDM_CD -- dragging in \c DM_Contract, \c DM_ContractBlocks, \c DM_RhoAtPoints
//! and the four HF J/K accumulators, none of which mixing has any business knowing about.  This is the
//! sibling of the \c GField / \c tFieldMixer split on the G-space side: each mixer family now names the
//! narrowest subject that supports its arithmetic.
//!
//! NB the two operations are BINARY, so they dispatch on both operands -- which C++ cannot do (single
//! dispatch only; the classic binary-method problem).  Every implementation therefore opens with a
//! same-kind \c dynamic_cast of \a other.  That is NOT a cost of this narrowing: the casts were already
//! there and are unavoidable (IrrepCD needs the other IrrepCD's MATRIX), and abstract-base to
//! abstract-base is exactly the cast this project sanctions.  A double-dispatch seam would remove them;
//! parked as not worth the second virtual hop today.
template <class T> class tMixableDensity
: public virtual tChargeDensity<T>
{
public:
    virtual void   MixIn        (const tMixableDensity<T>&,double)      =0;  //!< this = (1-c)*this + c*that
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const=0;  //!< convergence check
};

//  A charge density represented BY A DENSITY MATRIX: adds the matrix-only capabilities -- operator
//  contraction (DM_Contract) and the Hartree-Fock J/K accumulators.  Mixing is NOT here: it needs no
//  matrix, so it lives one level up on tMixableDensity (ISP).  Densities with no matrix (fits) are
//  tChargeDensity instead.
//
template <class T> class tDM_CD
: public virtual tMixableDensity<T>
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

    //! \brief \f$\rho\f$ at MANY points from caller-supplied basis tables -- the FIELD dual of
    //! DM_ContractBlocks.  \a Phi maps BasisSetID -> the (npoints x n) table \f$\Phi_{gi}=\chi_i(r_g)\f$
    //! for that block; the return is \f$\rho(r_g)=\sum_b\mathrm{Re}\,[\Phi_b D_b\Phi_b^\dagger]_{gg}\f$.
    //! The caller owns the (geometry-fixed, cacheable-for-the-whole-run) basis tables; the density keeps
    //! \f$D\f$ private and contracts them as GEMMs -- the seam that makes a mesh XC quadrature (Delta_XC)
    //! O(GEMM) per SCF iteration instead of re-evaluating Bloch sums pointwise.  A block whose ID is absent
    //! from \a Phi self-evaluates pointwise (correct, slower -- heals a caller's first pass before it has
    //! met every block).  Default: the plain pointwise loop over this density's ScalarFunction face, so
    //! EVERY density answers correctly; composite/leaf overrides accelerate.
    virtual rvec_t DM_RhoAtPoints(const rvec3vec_t& r, const std::map<std::string,mat_t<T>>& /*Phi*/) const
    {
        rvec_t ro(r.size());
        for (size_t g=0; g<r.size(); g++) ro[g]=(*this)(r[g]);
        return ro;
    }

    // HF/fit-specific (real, Gaussian-basis) -- the dcmplx plane-wave density NA-asserts these.
    // (The single-block diagonal contraction is NOT here: it is a private IrrepCD detail the whole-system
    //  AccumulateDirectAll loop reaches via AccumulateDirectBoth's self-pair branch -- see IrrepCD.)

    //! Whole-system HF Coulomb exploiting the ERI4 bra-ket symmetry \f$J(i,j)=J(j,i)^\mathsf{T}\f$
    //! (doc/ERI4Rework.md \S4/\S5.4): one uniform loop over canonical ab-irrep pairs \f$k\le l\f$ feeds BOTH
    //! Fock blocks (diagonal included -- see AccumulateDirectBoth), so J(j,i) is never built or stored.
    //! \a Jall is one block per irrep (same order as this composite's leaves) and pre-zeroed by the caller
    //! (the Vee term).  Only a composite / polarized density -- which spans every irrep block -- implements
    //! this; a lone leaf and the periodic path assert out.
    virtual void AccumulateDirectAll(std::vector<hmat_t<T>>& Jall) const
    { assert(false && "AccumulateDirectAll: only a composite/polarized density spans all irrep blocks"); }
    //! Contract ONE canonical pair (this=irrep i, \a other=irrep j, i<=j) into their Fock block(s) -- the
    //! SOLE per-pair helper the composite's AccumulateDirectAll loop calls.  Off-diagonal: fused scatter into
    //! both Ji and Jj.  Diagonal (other IS this): a single localized contraction into Ji (no bra-ket
    //! partner).  Only a leaf (IrrepCD) is a pair partner; composite/polarized inherit this asserting default.
    virtual void AccumulateDirectBoth(hmat_t<T>& Ji, hmat_t<T>& Jj, const tDM_CD<T>& other) const
    { assert(false && "AccumulateDirectBoth: only a leaf (irrep) density is a bra-ket pair partner"); }
    //! Exchange analogues of the two above.  Exchange is SAME-SPIN, so the whole-system build is per single-
    //! spin density (Polarized_CD sums the two spin channels' K into the same blocks for the RHF term).
    virtual void AccumulateExchangeAll(std::vector<hmat_t<T>>& Kall) const
    { assert(false && "AccumulateExchangeAll: only a composite/polarized density spans all irrep blocks"); }
    virtual void AccumulateExchangeBoth(hmat_t<T>& Ki, hmat_t<T>& Kj, const tDM_CD<T>& other) const
    { assert(false && "AccumulateExchangeBoth: only a leaf (irrep) density is a bra-ket pair partner"); }
};

using rChargeDensity = tChargeDensity<double>;  using cChargeDensity = tChargeDensity<dcmplx>;
using rDM_CD = tDM_CD<double>;
using cDM_CD = tDM_CD<dcmplx>;

//---------------------------------------------------------------------------------------
//
//  Store spin up and spin down as a ChargeDensity
//  Generic: Could be fitted or exact.
//  Templated on the matrix element type T like tComposite_CD: the <double> alias preserves the
//  molecular callers; the <dcmplx> instantiation is the polarized plane-wave (Bloch) density
//  (SymmetryUpgradePlan §4 tier 4b).  The periodic face (FourierDensityBase) forwards each accessor
//  as the ↑+↓ sum, so the total-density consumers (Hartree) see one density; the spin-native XC
//  terms reach the channels through GetChargeDensity(Spin).
//
template <class T> class tPolarized_CD
    : public virtual tDM_CD<T>
    , public virtual tLineageTracked<T>        // Layer-2: this top-level density tracks its SCF lineage head
    , public virtual ProjectedDensityBase<T>   // finite/molecular: an AO-projectable density
    , public FourierDensityBase<T>             // periodic: G-space rho-tilde (the ↑+↓ sum); empty base for double
{
public:
    virtual       tDM_CD<T>* GetChargeDensity(const Spin&)      =0;
    virtual const tDM_CD<T>* GetChargeDensity(const Spin&) const=0;

    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;   // sum both spins

    virtual double GetTotalCharge() const;  // <ro>
    virtual double GetTotalSpin  () const;  // No UT coverage// <up>-<down>

    // The spin children are mutated together (MixIn/ReScale below touch both), so either child's serial
    // tracks the polarized density's freshness; forward to Up.
    virtual size_t Version() const {return GetChargeDensity(Spin::Up)->Version();}

    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const;
    virtual void AccumulateDirectAll  (std::vector<hmat_t<T>>& Jall) const;  // sum both spins (Coulomb)
    virtual void AccumulateExchangeAll(std::vector<hmat_t<T>>& Kall) const;  // sum both spins (RHF exchange)

    virtual void   ReScale      (double factor              )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const tMixableDensity<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const;  //Convergence check.

    virtual double operator()(const rvec3_t&) const; // No UT coverage
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

    // Periodic (dcmplx) face -- the ↑+↓ sums of the channels' G-space/raster views (declared for both T
    // like tComposite_CD; the double bodies NA-assert).  Each channel composite already star-averages, so
    // the sums stay IBZ-symmetrized.
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const;   // empty if either channel lacks it
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const;
};

using Polarized_CD  = tPolarized_CD<double>;   // the molecular alias (source-compatible)
using cPolarized_CD = tPolarized_CD<dcmplx>;   // the polarized plane-wave (Bloch) density

//---------------------------------------------------------------------------------------
//
//  Capability face: a COLLINEAR two-channel spin-resolved density WITHOUT the tDM_CD (matrix)
//  contract -- the matrix-free polarized sibling of tPolarized_CD's channel accessor.  A spin-SAD
//  seed (doc/SCFSeedingPlan.md §10) has per-channel densities but no density matrix, so it cannot
//  be a tPolarized_CD (whose channels are tDM_CD, with the DM-only pure virtuals); it exposes its
//  channels through THIS face instead, and a spin-native consumer (XC_GridEngine::RhoPol) cross-casts
//  abstract->abstract and reads each channel through the plain tChargeDensity face (EvalBatch) --
//  capabilities live only on the types that have them (no asserting DM stubs).
//
template <class T> class tSpinResolved_CD
{
public:
    virtual ~tSpinResolved_CD() {}
    virtual const tChargeDensity<T>* GetChannel(const Spin&) const=0;   //!< Up/Down only (no None)
};

using rSpinResolved_CD = tSpinResolved_CD<double>;
using cSpinResolved_CD = tSpinResolved_CD<dcmplx>;

//! The magnetization ρ↑−ρ↓ as a real ScalarFunction (any T -- ρ_σ(r) is real for both lineages).
template <class T> class tSpinDensity : public virtual ScalarFunction<double>
{
public:
    tSpinDensity(tDM_CD<T>* up,tDM_CD<T>* down);
    ~tSpinDensity();
    virtual double operator()(const rvec3_t&) const; // No UT coverage
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage
private:
    tDM_CD<T>* itsSpinUpCD;
    tDM_CD<T>* itsSpinDownCD;
};

using SpinDensity = tSpinDensity<double>;

} //namespace