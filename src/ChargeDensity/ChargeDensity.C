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
export import qchem.Symmetry.Irrep;   // Irrep: the block identity (Phi-table cache key, basis-side)
export import qchem.ChargeDensity.FourierDensity;   // FourierDensityBase<T> (tPolarized_CD's periodic face)
import qchem.ScalarFunction;
export import qchem.Fitting.FunctionFitter;   // Fitting::ScalarProjector -- named on tDM_CD's own face,
                                              // so every consumer of that face needs it
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

//======================================================================================================
//  EXACT-EXCHANGE (HF) FACES -- V1.6/V1.8
//
//  These four operations used to be asserting defaults on tDM_CD<T>.  They are REAL-ONLY BY CONSTRUCTION:
//  Vee/Vxc are added by the molecular HF Hamiltonians alone, while the periodic Ham_PW_DFT adds
//  Vee_Hartree and never exact exchange -- which is why IrrepCD<dcmplx> carried four overrides saying "HF
//  is not applicable to a complex plane-wave density".  So they move to faces the REAL path inherits and
//  the complex path does not: the dcmplx instantiation now grows NOTHING, instead of one more denial per
//  method (the R2.8 smell).
//
//  T-TYPED, not type-erased, per doc/RealComplexPlan.md: an impossible pairing must fail to COMPILE
//  rather than throw.  Both faces mirror tMixableDensity<T>'s shape for that reason.
//
//! \brief The WHOLE-SYSTEM exact-exchange face: a density that spans EVERY irrep block, so it can drive the
//! canonical-pair sweep.  A composite/polarized density has it; a lone leaf does not -- and now cannot be
//! asked, where before it silently answered with a zeroed J.
template <class T> class tHF_System_CD
{
public:
    virtual ~tHF_System_CD() = default;
    //! Whole-system HF Coulomb via the ERI4 bra-ket symmetry \f$J(i,j)=J(j,i)^\mathsf{T}\f$
    //! (doc/ERI4Rework.md \S4/\S5.4): one loop over canonical pairs \f$k\le l\f$ feeds BOTH Fock blocks,
    //! so \f$J(j,i)\f$ is never built or stored.  \a Jall is one block per irrep, pre-zeroed by the caller.
    virtual void AccumulateDirectAll  (std::vector<hmat_t<T>>& Jall) const=0;
    //! Exchange analogue.  Exchange is SAME-SPIN, so the whole-system build runs per single-spin density.
    virtual void AccumulateExchangeAll(std::vector<hmat_t<T>>& Kall) const=0;
};

//! \brief The PAIR-PARTNER face: the composite<->leaf collaboration protocol for ONE canonical irrep pair.
//! Only a leaf is a pair partner, and this is the whole of what the composite needs from one -- it was never
//! public-face business (its only callers are the two loops in Imp/CompositeCD.C).
template <class T> class tHF_Pair_CD
{
public:
    virtual ~tHF_Pair_CD() = default;
    //! Contract the canonical pair (this = irrep i, \a other = irrep j, i<=j) into their Fock blocks.
    //! Off-diagonal: one fused scatter feeds both.  Diagonal (\a other IS this): a single localized
    //! contraction into \a Ji, since there is no bra-ket partner.
    virtual void AccumulateDirectBoth  (hmat_t<T>& Ji, hmat_t<T>& Jj, const tHF_Pair_CD<T>& other) const=0;
    virtual void AccumulateExchangeBoth(hmat_t<T>& Ki, hmat_t<T>& Kj, const tHF_Pair_CD<T>& other) const=0;

    // The two below are called by a PEER partner, not by clients -- protected would not do, since a sibling
    // instance is not `this`.  Public on a face that is itself the composite<->leaf protocol is honest: the
    // face's whole existence is that collaboration, and nothing outside it holds a tHF_Pair_CD.
    //! \brief The partner's half, named for the OPERATION rather than the data (V1.8).
    //!
    //! The pair contraction belongs to the BASIS -- the caller ends up handing (D_i, bs_i, D_j, bs_j) to
    //! \c Orbital_HF_IBS::AccumulateDirectBoth -- so a partner's only job is to present ITS OWN block to
    //! that routine, given the caller's.  Expressed as "finish this contraction with your block", NOT as
    //! "give me your density matrix": \c IrrepCD deliberately has no \c GetDensityMatrix(), and CLAUDE.md
    //! cites that very class as the example of preferring a high-level operation to an exposed member.  An
    //! abstract block ACCESSOR would have smuggled the getter back in behind an interface.
    //!
    //! This is what retires V1.8's abstract->concrete `IrrepCD*` cast: the caller no longer needs the
    //! partner's TYPE, only its cooperation.
    virtual void CompleteDirectPair  (hmat_t<T>& Ji, hmat_t<T>& Jj, const hmat_t<T>& Di,
                                      const BasisSet::Orbital_HF_IBS<T>* bs_i) const=0;
    virtual void CompleteExchangePair(hmat_t<T>& Ki, hmat_t<T>& Kj, const hmat_t<T>& Di,
                                      const BasisSet::Orbital_HF_IBS<T>* bs_i) const=0;
};

//! Empty stand-ins: the periodic (dcmplx) path has no exact exchange, so it inherits nothing at all.
struct NoHF_System {};
struct NoHF_Pair   {};
//! The HF faces on the finite path (T=double), empty bases on the periodic path (T=dcmplx) -- the same
//! idiom as ProjectedDensityBase / FourierDensityBase directly below.
// (No HF_SystemBase alias: the two whole-system densities carry their sweep in a CRTP mixin of their own
// -- Composite_HFSystem / Polarized_HFSystem -- so the face is inherited by those, real path only.)
template <class T> using HF_PairBase   = std::conditional_t<std::is_same_v<T,double>, tHF_Pair_CD<double>,   NoHF_Pair>;

//! The COULOMB-metric projection face for the finite path (T=double), the empty base for the periodic path
//! (T=dcmplx).  V1.16: a matrix-backed density HAS the Coulomb RHS, so it takes the face that requires it --
//! rather than inheriting a metric DEFAULT plus a poisoned sibling it had to remember not to call.  The
//! matrix-FREE seeds take Fitting::OverlapMetric_ProjectedDensity instead (see NumericCD/SeedCD).
template <class T> using ProjectedDensityBase =
    std::conditional_t<std::is_same_v<T,double>, Fitting::CoulombMetric_ProjectedDensity, NoProjectedDensity>;
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
    //! \brief The ENERGY matrix of a density-dependent term: the block whose D-contraction IS this term's
    //! energy contribution, \f$E=\mathrm{Tr}(D\,B)\f$.  \a cd is the density defining the term (the DFT-only
    //! face -- a fit has no matrix); the basis/spin identify which irrep block the caller wants.
    //!
    //! DELIBERATELY NOT SPELLED \c GetMatrix.  For nearly every term the energy matrix IS the Fock block
    //! (\f$E=D\cdot V\f$), and \c tDynamic_HT's default forwards to it.  The xc family is the exception --
    //! \f$v_{xc}=\delta E_{xc}/\delta\rho=\epsilon_{xc}+\rho\,\partial\epsilon_{xc}/\partial\rho\f$, so
    //! \f$\int v_{xc}\rho\neq\int\epsilon_{xc}\rho=E_{xc}\f$ (Slater exchange: a factor 4/3) -- and it
    //! overrides this with the ENERGY-DENSITY block \f$\langle i|\tilde\epsilon_{xc}|j\rangle\f$.  While the
    //! two faces shared the \c GetMatrix spelling a term could only ever offer ONE of the two matrices, which
    //! is why the xc terms had to hang a second \c tDynamic_CC object (an "eps adapter") off the side purely
    //! to carry the other one.  Distinct names dissolve that collision (V1.3).
    virtual const hmat_t<T>& GetEMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*) const=0;
};

//====================================================================================================
//  REAL-BLOCK CONTRACT CLIENTS (doc/RealComplexPlan.md Step 3c-2b).  A real TRIM block inside a
//  COMPLEX run energy-contracts against the run's terms.  The STATIC case needs no new face at all:
//  Static_HT_RealBlock::GetMatrix has exactly tStatic_CC<double>'s signature, so the term's real-block
//  mixin simply IS a tStatic_CC<double> and the real leaf's native DM_Contract serves.  The DYNAMIC
//  case cannot reuse tDynamic_CC<double> -- its energy matrix takes the RUN's density, which is
//  complex-faced -- hence this one run-typed client, implemented by the dynamic term's real-block
//  mixin (default: the E=D·V identity, i.e. its real GetMatrix; the periodic xc terms never take the
//  DM_Contract route, so no eps override is needed).
class Dynamic_CC_RealBlock
{
public:
    virtual ~Dynamic_CC_RealBlock() {};
    virtual const hmat_t<double>& GetEMatrixR(const tobs_t<double>*,const Spin&,
                                              const tChargeDensity<dcmplx>*) const=0;
};

//! \brief The capability of a REAL block density living inside a COMPLEX run: energy-contract this
//! block's (real) D against the run-typed dynamic client above.  Carried by the real periodic leaf
//! alone (PeriodicIrrepCD<double>); the composite's cross-scalar contract arm cross-casts to it.
class RealBlockEnergy_CD
{
public:
    virtual ~RealBlockEnergy_CD() {};
    virtual double DM_ContractE(const Dynamic_CC_RealBlock*, const tChargeDensity<dcmplx>*) const=0;
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

    //  rho at MANY points is INHERITED, not re-declared: ScalarFunction<double>::operator()(rvec3vec_t)
    //  already IS that function (same signature after alias expansion), with the same pointwise-loop
    //  default.  A density that can vectorize the batch overrides THAT (FourierMixCD's factorized phase
    //  tables); the DM densities have the richer, table-driven ProjectOnto instead.  The mesh-XC
    //  quadrature engine's entry for non-DM densities.  (There WAS a duplicate `EvalBatch` spelling here,
    //  and it had forked: FourierMixCD overrode only EvalBatch, so any caller reaching a density through
    //  the neutral ScalarFunction batch face would have silently taken the per-G std::exp loop instead.
    //  Every such caller happened to spell it EvalBatch -- the fast path was correct by luck of spelling.
    //  R1.5.)

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
//! yet it used to take a \c tDM_CD -- dragging in \c DM_Contract, \c DM_ContractBlocks, \c ProjectOnto
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
    //! cached blocks WITHOUT a per-irrep GetMatrix round-trip through DM_Contract.  Every density
    //! implements it, including the periodic path (PW_Hartree contracts its long-range blocks here).
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const=0;

    //! \brief PROJECT ME ONTO \a p's fit basis -- the FIELD dual of \c DM_ContractBlocks:
    //! \f$c_a=\langle f_a|\rho\rangle/\langle f_a|f_a\rangle\f$, realized here as
    //! \f$\sum_b\mathrm{Re}\,[\Phi_b D_b\Phi_b^\dagger]_{gg}\f$.
    //!
    //! THE DENSITY ASKS THE FITTER (2026-08-24).  \f$D\f$ is this class's private business and the
    //! \f$\Phi\f$ TABLE is the fit basis's, so neither travels: the fitter hands over the CONTRACTIBLE
    //! integral (\c Projector3, a shallow handle on its cached table) and each block contracts its own
    //! \f$D\f$ or thin \f$L\f$ into it.  MIXED-aware for free -- a real TRIM block asks with \c double
    //! while its complex siblings ask with \c dcmplx (doc/RealComplexPlan.md 3c-3).
    //!
    //! \note WHY NOT \c operator()(rvec3vec_t) (user, 2026-08-25).  That face receives COORDINATES, and
    //! the table \f$\chi_i(r_g)\f$ is not reachable from them -- so an override there could only
    //! RECOMPUTE it, which costs what the pointwise sweep costs (measured: 80 ms vs 77 ms at 4000×16; the
    //! ~500× is the CACHE, not the contraction).  Nor could it cache: a density is a fresh object every
    //! SCF iteration (\c TOrbitalsImp::GetChargeDensity news one), and it is asked for \f$\rho\f$ once per
    //! iteration, so a density-side table would have a zero hit rate.  \a p is exactly that face PLUS the
    //! identity needed to reach the run-lifetime cache.  Renamed from \c DM_RhoAtPoints for the same
    //! reason: it does not return values at points, and "points" was the last of the vocabulary the
    //! 2026-08-23 pass took off everything else.
    //!
    //! \return MY EXPANSION COEFFICIENTS over \a p's fit basis -- one per fit FUNCTION, not one per point.  For a
    //! \f$\delta\f$ basis the two coincide numerically (\f$c_g=\rho(r_g)\f$), which is exactly why the
    //! pointwise-nonlinear XC functional may be applied to them directly; the name of the array is the
    //! coefficient vector, and the equality is a property of the representation.
    //!
    //! PURE, with no default (2026-08-25).  There WAS one -- hand \a p this density's own
    //! \c ScalarFunction face and let the fitter sample it -- and it was DEAD: measured 0 calls across
    //! all 760 tests, with a control probe on the \c IrrepCD override firing to prove the instrument
    //! worked.  It could not be otherwise: this face is \c tDM_CD, so reaching it means HAVING a matrix,
    //! and a density that has one always has a better answer than being sampled pointwise.  A density with
    //! NO matrix never arrives here at all -- it is not a \c tDM_CD, and its caller asks
    //! \c ScalarProjector::Project directly.  So the field route was a default on the wrong face.
    //!
    //! Same shape as V1.16 on the Coulomb side: each capability is a face you either HAVE or do not, and
    //! there is no fallback left to hit by accident.
    virtual rvec_t ProjectOnto(const Fitting::ScalarProjector& p) const=0;

    // The exact-exchange (HF) accumulators are NOT here any more -- see tHF_System_CD / tHF_Pair_CD
    // below (V1.6).  They were four asserting defaults on this general face, which every concrete family
    // relied on for exactly two of the four, and which the dcmplx instantiation overrode with four more
    // asserts saying "HF does not apply to me".  Two costs: the void assert-only bodies are silent NO-OPS
    // under -DNDEBUG, so a bare leaf reaching Vee::AccumulateAll yielded a ZEROED J and a silently wrong
    // Fock in the build we actually ship; and every future density kind had to re-declare the same denials.

    //! \brief The WHOLE-SYSTEM Fock route for THIS block's basis, or null if it has none (V1.31).
    //!
    //! A capability PROBE in the \c isPeriodicCell sense (R2.4): the caller branches on presence, ONCE, at
    //! the top of the sweep.  Non-null means the basis builds one Fock from the total AO density and slices
    //! it per irrep, so the canonical-PAIR loop above is the wrong shape for it -- a SALC basis has no
    //! per-irrep-pair ERI4 blocks at all (R1.7), and driving the pair loop made it rebuild the same
    //! whole-AO Fock ~N times.  Null (every ERI4 basis) means the pair loop applies, unchanged.
    virtual const BasisSet::WholeSystemFock_IBS<T>* WholeSystemFock() const {return nullptr;}
    //! Fold THIS block's density into the whole-system (AO) total.  Only meaningful on a leaf whose basis
    //! answered \c WholeSystemFock(); the composite calls it once per leaf before the single build.
    virtual void AddAODensity(hmat_t<T>& Dao) const
    { assert(false && "AddAODensity: only a leaf (irrep) density carries a block density matrix"); }
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
//! \brief The reciprocal-space trio for a POLARIZED density -- periodic path only (V1.7): the ↑+↓ sums of
//! the two channels' G-space/raster views.  Each channel composite already star-averages, so the sums stay
//! IBZ-symmetrized.  CRTP, like its leaf and composite siblings.
template <class Pol> class Polarized_Fourier : public virtual FourierDensity
{
public:
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const;   // empty if either channel lacks it
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const;
private:
    const Pol& self() const {return static_cast<const Pol&>(*this);}
};

template <class T, class Pol> using PolarizedFourierBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, Polarized_Fourier<Pol>, NoFourierDensity>;

//! \brief The whole-system exact-exchange sweep for a POLARIZED density -- real path only (V1.6,
//! completing it): each channel spans every block, so both simply drive their own sweep into the shared
//! Fock blocks.  CRTP like its composite and leaf siblings, and for the same reason: inheriting the
//! \c tHF_System_CD face conditionally while DECLARING its methods unconditionally left the dcmplx
//! instantiation carrying members that override nothing and can only throw.
template <class Pol> class Polarized_HFSystem : public virtual tHF_System_CD<double>
{
public:
    virtual void AccumulateDirectAll  (std::vector<hmat_t<double>>& Jall) const;
    virtual void AccumulateExchangeAll(std::vector<hmat_t<double>>& Kall) const;
private:
    const Pol& self() const {return static_cast<const Pol&>(*this);}
};

template <class T, class Pol> using PolarizedHFBase =
    std::conditional_t<std::is_same_v<T,double>, Polarized_HFSystem<Pol>, NoHF_System>;

template <class T> class tPolarized_CD
    : public virtual tDM_CD<T>
    , public virtual tLineageTracked<T>        // Layer-2: this top-level density tracks its SCF lineage head
    , public virtual ProjectedDensityBase<T>   // finite/molecular: an AO-projectable density
    , public PolarizedFourierBase<T,tPolarized_CD<T>>   // reciprocal trio: periodic path only (V1.7)
    , public PolarizedHFBase<T,tPolarized_CD<T>>   // whole-system exact exchange: real path only (V1.6)
{
public:
    virtual       tDM_CD<T>* GetChargeDensity(const Spin&)      =0;
    virtual const tDM_CD<T>* GetChargeDensity(const Spin&) const=0;

    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;   // sum both spins
    //! \copydoc tDM_CD::ProjectOnto
    //! BOTH CHANNELS SUMMED -- \f$\rho=\rho_\uparrow+\rho_\downarrow\f$, the same shape as
    //! \c DM_ContractBlocks above and as \c Polarized_Fourier::GetRhoOnGrid on the raster side.
    //! Stated rather than inherited since 2026-08-25: the base default was measured dead and removed, and
    //! a polarized density genuinely does have an answer here (a spin-agnostic consumer asking a polarized
    //! density for \f$\rho\f$ wants the total).  Note the SPIN-RESOLVED consumer does not come here at
    //! all -- \c XC_SinglesQuadrature::RhoPol asks each CHANNEL, because it needs them apart.
    virtual rvec_t ProjectOnto(const Fitting::ScalarProjector&) const;

    virtual double GetTotalCharge() const;  // <ro>
    virtual double GetTotalSpin  () const;  // No UT coverage// <up>-<down>

    // The spin children are mutated together (MixIn/ReScale below touch both), so either child's serial
    // tracks the polarized density's freshness; forward to Up.
    virtual size_t Version() const {return GetChargeDensity(Spin::Up)->Version();}

    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const;
    // The whole-system J/K sweep is NOT declared here (V1.6 ISP): it lives in Polarized_HFSystem, which
    // only the REAL instantiation inherits -- the periodic path has no exact exchange to deny.

    virtual void   ReScale      (double factor              )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const tMixableDensity<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const;  //Convergence check.

    virtual double operator()(const rvec3_t&) const; // No UT coverage
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage

    // The reciprocal trio is NOT declared here (V1.7 ISP) -- see Polarized_Fourier, inherited on the
    // periodic path only.  It used to be declared for both T with the double bodies NA-asserting.
};

using Polarized_CD  = tPolarized_CD<double>;   // the molecular alias (source-compatible)
using cPolarized_CD = tPolarized_CD<dcmplx>;   // the polarized plane-wave (Bloch) density

//---------------------------------------------------------------------------------------
//
//  Capability face: a COLLINEAR two-channel spin-resolved density WITHOUT the tDM_CD (matrix)
//  contract -- the matrix-free polarized sibling of tPolarized_CD's channel accessor.  A spin-SAD
//  seed (doc/SCFSeedingPlan.md §10) has per-channel densities but no density matrix, so it cannot
//  be a tPolarized_CD (whose channels are tDM_CD, with the DM-only pure virtuals); it exposes its
//  channels through THIS face instead, and a spin-native consumer (XC_Quadrature::RhoPol) cross-casts
//  abstract->abstract and reads each channel through the plain tChargeDensity face (the batched op()) --
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

//---------------------------------------------------------------------------------------
//
//  Capability face: a FIELD-backed density that RETAINS THE DENSITY MATRIX IT WAS BUILT FROM.
//
//  WHY IT EXISTS.  Under ρ̃-mixing (Kerker/Pulay) the density driving the Fock build is a G-space FIELD
//  with no D, so a consumer that wants ρ(r) must inverse-transform a TRUNCATED Fourier series at every
//  point.  For the Hartree term that is exactly right -- Poisson is LINEAR and diagonal in G, so the
//  preconditioned field is the object the operator wants.  For XC it is wrong twice over, and both were
//  MEASURED on the MnO benchmark (2026-08-21):
//    * COST: the direct summation is O(npts x nG) with no BLAS -- 5.19 s per sampling, 51% of the whole
//      run -- against 0.042 s for the same rho through the (factored) DM GEMM.
//    * ACCURACY: a truncated series RINGS, and it rings NEGATIVE where the true density is sharpest, at
//      the nuclear cusps the atom-centred mesh exists to integrate.  8.5-18% of mesh points came back
//      with rho<0, min rho ~ -0.14 -- and the functionals guard `if (rho>0)`, so every one of those
//      points contributes ZERO to v_xc and E_xc with no diagnostic.  rho = ||L^dag Phi||^2 cannot.
//  Hartree and XC therefore want DIFFERENT REPRESENTATIONS OF THE SAME DENSITY, and this face is how the
//  second one stays reachable: the mixer retains the DM-backed density it built the field from, and a
//  quadrature consumer cross-casts for it.  At the fixed point the two agree, so this changes the SCF
//  TRAJECTORY, not the answer.
//
//  SHARED ownership, deliberately: the retained density must outlive the Mix() call that produced it (XC
//  samples it later in the iteration), which is exactly the case tDensityMixer's plain-reference subject
//  does NOT cover.  A raw back-pointer here would be valid only between mixes and its failure mode --
//  a stale D silently producing a wrong V_xc -- is invisible in Release.
//
template <class T> class tDM_Sourced_CD
{
public:
    virtual ~tDM_Sourced_CD() {}
    //! The DM-backed density this field was built from; null when none has been deposited (iteration 0,
    //! or a mixer that does not track one) -- callers fall back to the field's own \c operator().
    virtual std::shared_ptr<const tDM_CD<T>> DMSource() const=0;

    //! \brief THE EFFECTIVE MIXING FRACTION this field actually applied -- the FRACTION OF THE UPDATE that
    //! survived the preconditioner, \f$\alpha_{\rm eff}=\lVert\tilde\rho_{mix}-\tilde\rho_{in}\rVert_2 /
    //! \lVert\tilde\rho_{out}-\tilde\rho_{in}\rVert_2\f$.
    //!
    //! WHY IT BELONGS ON THIS FACE.  A consumer that reaches around the field to \c DMSource() gets the
    //! UNDAMPED output density, so it has stepped outside the mixing the rest of the SCF is subject to.
    //! Feeding XC that raw while Hartree keeps the preconditioned field is a HALF-DAMPED map, and it
    //! measurably wrecks convergence (NaF DIIS 34 -> 100+ iterations, 2026-08-21).  So the two pieces are
    //! wanted TOGETHER, always: the exact density, and how much of it to take.
    //!
    //! MEASURED, NOT MODELLED.  Being a ratio of norms of objects the mixer already holds, it needs no
    //! knowledge of WHICH preconditioner ran: for Kerker it evaluates to
    //! \f$\alpha\sqrt{\sum f^2|\delta\tilde\rho|^2/\sum|\delta\tilde\rho|^2}\f$ (α times the
    //! residual-power-weighted RMS of the filter, so it rises toward α as the residual moves to high G);
    //! for a plain linear mix it is exactly α; and for a Pulay step it may EXCEED 1, which is honest --
    //! extrapolation overshoots on purpose.
    //!
    //! \return 0 when nothing has been mixed yet (iteration 0) -- callers should treat that as "no damping
    //! information", not as "damp to zero".
    virtual double EffectiveAlpha() const=0;
};

using rDM_Sourced_CD = tDM_Sourced_CD<double>;
using cDM_Sourced_CD = tDM_Sourced_CD<dcmplx>;

//! The magnetization ρ↑−ρ↓ as a real ScalarFunction (any T -- ρ_σ(r) is real for both lineages).
template <class T> class tSpinDensity : public virtual ScalarFunction<double>
{
public:
    //! TAKES OWNERSHIP of both channels (V1.25).  The old form held raw pointers AND deleted them in a
    //! hand-written dtor while declaring no copy ctor or assignment -- so copying it was a double-delete.
    //! It was safe only because its one caller happened to hand it two freshly-built densities.
    tSpinDensity(std::unique_ptr<tDM_CD<T>> up, std::unique_ptr<tDM_CD<T>> down);
    virtual double operator()(const rvec3_t&) const; // No UT coverage
    virtual rvec3_t  Gradient  (const rvec3_t&) const; // No UT coverage
private:
    std::unique_ptr<tDM_CD<T>> itsSpinUpCD;
    std::unique_ptr<tDM_CD<T>> itsSpinDownCD;
};

using SpinDensity = tSpinDensity<double>;

} //namespace