// File: CompositeCD.C  Composite charge density: ONE set of irrep block densities over the FULL Irrep (spatial ⊗ Spin).
module;
#include <vector>
#include <memory>
#include <cstddef>
#include <map>
#include <string>
#include <variant>
export module qchem.CompositeCD;
export import qchem.ChargeDensity;
export import qchem.ChargeDensity.FourierDensity;   // G-space rho-tilde (summed over k-blocks)
export import qchem.Symmetry.Lattice_3D.SpaceGroup; // ReciprocalOp {U|τ} -- the IBZ density-symmetrization ops (glide phase)
export import qchem.Symmetry.Irrep;                 // the block label (spatial ⊗ ms)
import qchem.ChargeDensity.Types;


export namespace qchem::ChargeDensity
{

//! \brief THE CHILD SLOT (doc/RealComplexPlan.md §4, Step 2): one block density, typed by ITS OWN scalar
//! rather than the composite's face -- so children of one composite may differ (a real TRIM block beside
//! general-k complex blocks on a mixed mesh).  The alternatives are the ABSTRACT per-block face, ownership
//! included.  Aggregation stays single-source: scalar-independent operations visit with ONE generic lambda;
//! the T-typed operations forward to the same-scalar alternative (the cross arm becomes reachable -- and
//! gets its narrowing implementation -- when Step 3 un-pins the basis type per block).
using cd_child_t = std::variant<std::unique_ptr<tDM_CD<double>>, std::unique_ptr<tDM_CD<dcmplx>>>;
//! Non-owning mirror of the child slot: what the walk order and the channel VIEWS hold (V1.37).
using cd_ref_t   = std::variant<tDM_CD<double>*, tDM_CD<dcmplx>*>;

//--------------------------------------------------------------------------
//
//  Full charge density represented Compositely as sum of density matricies.
//  Templated on the matrix element type T (rX/cX); the <double> alias preserves the existing
//  real callers, the <dcmplx> instantiation aggregates the plane-wave (Bloch-irrep) densities.
//
//! \brief The reciprocal-space trio for a COMPOSITE -- periodic path only (V1.7).  Each block already
//! carries its BZ weight, so every one of these is a straight sum over the leaves.  CRTP for the same reason
//! as the leaf's: reach the contained blocks without exposing them.
template <class Comp> class Composite_Fourier : public virtual FourierDensity
{
public:
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const;   // empty if any block lacks the raw path
    virtual ΔG_Map GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const;
    //! The raw/average pair (see FourierDensity): this composite owns the crystal point ops, so it owns
    //! the star average -- split out so the polarized case merges BOTH channels raw and averages ONCE.
    virtual ΔG_Map GetRepulsion3C_Raw(const BasisSet::cFIT_CD_ABS& c) const;
    virtual void   StarAverage(ΔG_Map& rg) const;
private:
    const Comp& self() const {return static_cast<const Comp&>(*this);}
};

template <class T, class Comp> using CompositeFourierBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, Composite_Fourier<Comp>, NoFourierDensity>;

//! \brief The whole-system exact-exchange sweep for a COMPOSITE -- real path only (V1.6, completing it).
//! Same CRTP shape as \c Composite_Fourier, and for the same reason the LEAF got one: inheriting the
//! \c tHF_System_CD face conditionally but DECLARING its two methods unconditionally left
//! \c tComposite_CD<dcmplx> carrying exact-exchange members that override nothing and can only throw --
//! the very "denial the other half of the hierarchy must write" this item set out to remove, surviving
//! one level above the leaf.  Declared here, the complex composite grows nothing at all.
template <class Comp> class Composite_HFSystem : public virtual tHF_System_CD<double>
{
public:
    virtual void AccumulateDirectAll  (std::vector<hmat_t<double>>& Jall) const;
    virtual void AccumulateExchangeAll(std::vector<hmat_t<double>>& Kall) const;
private:
    const Comp& self() const {return static_cast<const Comp&>(*this);}
    //! One spin group's sweep (blocks begin..end of the walk) into \a Fall -- Coulomb or exchange.
    void SweepGroup(size_t begin, size_t end, std::vector<hmat_t<double>>& Fall, bool exchange) const;
};

template <class T, class Comp> using CompositeHFBase =
    std::conditional_t<std::is_same_v<T,double>, Composite_HFSystem<Comp>, NoHF_System>;

//! \brief ONE composite over FULL Irreps (spatial ⊗ Spin) -- the SCF density of every run, polarized or not
//! (doc/CleanupCandidates.md V1.37).
//!
//! THE RULING.  Pol/UnPol is the IMPOSED SPIN SUBGROUP (\c qchem::SpinGroup), the same kind of decision as
//! imposing a point group -- so it is a property of the LABELS, never of the container.  A polarized run
//! inserts blocks labelled \c Spin::Up and \c Spin::Down; an unpolarized run inserts \c Spin::None blocks,
//! each of which is the FOLDED DOUBLET and counts for two because its irrep says so
//! (\c Irrep::GetDegeneracy; the fold rides in the block's occupations, so this class applies no weight).
//! The two-level \c Polarized{Composite,Composite} tree this replaced hard-coded two collinear channels into
//! a type; a flat set over double-group irreps (spin INSIDE G) needs no new container.
//!
//! THE SPIN STRUCTURE IS A VIEW.  \c GetChannel(s) answers the \c tSpinResolved_CD face with a composite over
//! the blocks whose irrep carries \a s -- non-owning, cached, alive as long as this density -- and carries
//! every face this one does (matrix, Fourier, HF sweep), so a spin-native consumer (XC, exchange, the
//! per-channel mixer) reads a channel exactly as a spin-agnostic one (Hartree, 1E) reads the total.  A
//! composite that resolves no such spin answers null (the R2.4 probe idiom).
//!
//! THE WALK ORDER IS THE BASIS ORDER.  Blocks are kept in insertion order -- the order \c tCompositeWF
//! builds them, i.e. \c itsBS->Iterate per spin irrep -- because the whole-system HF sweep pairs block k
//! of a SPIN GROUP with Fock block \c Jall[k] positionally (exchange is same-spin; the Fock blocks are per
//! spatial irrep).  Every other aggregation is a plain sum over all blocks.  (V1.37 first landed with the
//! sums grouped by spin to reproduce the old ↑-sum + ↓-sum tree bit for bit; user ruling 2026-09-14: clean
//! code over a 1e-16 reordering -- the flat sum stayed, totals unchanged at printed precision, 857/857.)
template <class T> class tComposite_CD
    : public virtual tDM_CD<T>
    , public virtual tLineageTracked<T> // Layer-2: this top-level density tracks its SCF lineage head
    , public virtual tSpinResolved_CD<T> // the channel VIEW face (V1.37)
    , public ProjectedDensityBase<T> // AO projection on the finite (double) path; empty on the periodic path
    , public CompositeFourierBase<T,tComposite_CD<T>>   // reciprocal trio: periodic path only (V1.7)
    , public CompositeHFBase<T,tComposite_CD<T>>        // whole-system exact exchange: real path only (V1.6)
{
public:
    //! \a pointOps = the reciprocal point group for IBZ density symmetrization (doc/GPWPlan1.md item 3).  The
    //! G-space density accessors then return the STAR AVERAGE, which is what makes an IBZ-reduced density exact
    //! (the star weights alone give only the correct band sum).  Default {} = trivial group {E} = exact no-op
    //! -- molecules / Γ / unreduced crystals pass through untouched (the general form; "no symmetry" = trivial).
    //! It is a ctor argument, not a setter: the symmetry is a fixed property of the density, set once at build.
    explicit tComposite_CD(std::vector<Symmetry::Lattice_3D::ReciprocalOp> pointOps = {});
    ~tComposite_CD();
    //! TAKES OWNERSHIP of \a cd (V1.25), labelled by its FULL irrep \a qns (spatial ⊗ ms -- the spin is
    //! what the channel views filter on).  TWO overloads, one per child scalar (doc/RealComplexPlan.md
    //! Step 2): the child slot is typed by the BLOCK, not by this composite's face, so either alternative
    //! may be inserted regardless of T.  A molecular composite simply never receives the complex one.
    void Insert(std::unique_ptr<tDM_CD<double>> cd, const Irrep& qns);
    void Insert(std::unique_ptr<tDM_CD<dcmplx>> cd, const Irrep& qns);

    //! tSpinResolved_CD: the \a s channel as a composite VIEW over this density's \a s blocks; null when no
    //! block carries \a s.  A single-spin composite (a view, or an unpolarized run's density asked for
    //! \c None's own channel) IS its channel and answers itself.
    virtual const tChargeDensity<T>* GetChannel(const Spin& s) const;

    // The whole-system J/K sweep is NOT declared here (V1.6 ISP): it lives in Composite_HFSystem, which
    // only the REAL instantiation inherits -- so tComposite_CD<dcmplx> declares nothing and defines nothing
    // for an operation the periodic path never has (Ham_PW_DFT adds Vee_Hartree, never exact exchange).
    virtual double DM_Contract(const tStatic_CC<T>*) const;
    virtual double DM_Contract(const tDynamic_CC<T>*,const tDM_CD<T>*) const;
    virtual double DM_ContractBlocks(const std::map<std::string,hmat_t<T>>&) const;   // sum over irrep blocks
    //! \copydoc tDM_CD::ProjectOnto
    //! ALL BLOCKS SUMMED -- for a polarized composite \f$\rho=\rho_\uparrow+\rho_\downarrow\f$, the total a
    //! spin-agnostic consumer wants.  The SPIN-RESOLVED consumer does not come here: it asks each CHANNEL.
    virtual rvec_t ProjectOnto(const Fitting::ScalarProjector&) const;

    virtual double GetTotalCharge      (                     ) const;

    // The blocks are mutated together (MixIn/ReScale fan out to all), so any block's serial tracks the
    // composite's freshness; forward to the first.  Empty composite -> 0 (the "no density yet" sentinel).
    // NB the Up channel view shares its first block with the total, hence its serial: a cache keyed on
    // Version() must never be warmed with the total and then asked for a channel (FittedVxcPol's note).
    virtual size_t Version() const
    {return itsBlocks.empty() ? 0 : std::visit([](const auto* c){return c->Version();}, itsBlocks.front().cd);}

    virtual double FitGetConstraint() const {return GetTotalCharge();}   // AO fit RHS: the charge N
    virtual rvec_t GetRepulsion3C(const BasisSet::rFIT_CD_ABS*) const;

    virtual void   ReScale      (double factor         )      ;  // No UT coverage//Ro *= factor
    virtual void   MixIn        (const tMixableDensity<T>&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const tMixableDensity<T>&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const rvec3_t&) const;
    virtual rvec3_t  Gradient  (const rvec3_t&) const;

    // The reciprocal trio is NOT declared here (V1.7 ISP) -- see Composite_Fourier below, inherited only
    // on the periodic path, so the finite composite is not asked questions it cannot answer.

private:
    friend class Composite_Fourier<tComposite_CD<T>>;   // sums the blocks; nothing is exposed publicly
    friend class Composite_HFSystem<tComposite_CD<T>>;  // drives the canonical-pair sweep over the blocks
    tComposite_CD(const tComposite_CD&);

    //! One block: its full label and the (non-owning) block density.
    struct Block { Irrep qns; cd_ref_t cd; };
    //! A VIEW: non-owning, over \a blocks, sharing the parent's point ops.  Only GetChannel builds one.
    tComposite_CD(std::vector<Block> blocks, const std::vector<Symmetry::Lattice_3D::ReciprocalOp>& pointOps);
    //! One contiguous run of same-spin blocks in the walk (the HF sweep's unit; also what the views filter on).
    struct Group { Spin s; size_t begin, end; };
    std::vector<Group> Groups() const;
    //! Rebuild the channel views after an Insert (eager, so const readers never race on a lazy build).
    void RebuildChannels();

    std::vector<cd_child_t> itsOwned;    // ownership (EMPTY in a view); the §4 child slot, per-block scalar
    std::vector<Block>      itsBlocks;   // the walk order: insertion == the basis order, per spin irrep
    std::vector<Symmetry::Lattice_3D::ReciprocalOp> itsPointOps;   //!< reciprocal {U|τ} point group for IBZ symmetrization ({E} when empty)
    std::map<Spin,std::unique_ptr<tComposite_CD>> itsChannels;     //!< the channel views (empty in a view / single-spin set)
};

using rComposite_CD = tComposite_CD<double>;   using cComposite_CD = tComposite_CD<dcmplx>;

} //namespace
