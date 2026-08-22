// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
module;
#include <cassert>   // AdoptMOMReference's OrbitalView cross-cast guard
#include <memory>
#include <functional>
#include <string>
#include <iosfwd>
#include <vector>
#include <optional>
export module qchem.SCFIterator;
import qchem.SCFIterator.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
import qchem.WaveFunction.SCF;
import qchem.ElectronConfiguration.OccupationPolicy;   // the iterator's occupation SLOT (V1.11 inc 3)
export import qchem.SCFParams;
export import qchem.ChargeDensity.Seed;   // SeedStrategy / MakeSeedDensity
import qchem.LASolver;   // qchem::Ortho (the basis-overlap orthogonalisation knob, forwarded to the WF)
import qchem.BasisSet.Fit_IBS;   // BasisSet::FIT_SF_ABS<T> (the G-space fit basis for Kerker rho-tilde extraction)
export import qchem.ChargeDensity.DensityMixer;   // tDensityMixer<T> (the density-face of SCF convergence)
import qchem.SCFIterator.LoopDriver;   // tLoopDriver<T> + Fixed/DirectMin concretes (the loop-face seam)

export using qchem::EnergyBreakdown;
using qchem::ChargeDensity::tDM_CD;
using qchem::ChargeDensity::tChargeDensity;

export namespace qchem::SCFIterator
{

//! Process-wide diagnostic toggle (default OFF), mirroring \c qchem::Hamiltonian::ReportGridCharge.  When
//! true, the verbose per-iteration SCF line appends the
//! frontier spectrum -- ε_HOMO, ε_LUMO and the gap ε_LUMO−ε_HOMO (Ha) -- from the current orbital energy
//! levels.  The band-gap instrument for the NaF Γ-instability (doc/GPWPlan §0b″): a near-degenerate
//! HOMO/LUMO (gap → 0) is the giant-response hypothesis; watching the gap per iteration lets the spurious
//! level be seen diving across the Fermi edge one step before each energy spike.  Flip it around
//! \c Iterate and reset it (it is a static, so it leaks between tests otherwise).
bool& ReportBandGap();

//! Per-iteration progress, handed to an Observer each SCF step so a client (a GUI,
//! a logger) can watch convergence live without owning the loop.  All real-valued
//! (energies/residuals are real for both the rX and cX paths).
struct SCFProgress
{
    size_t iteration;
    double energy;       //!< total energy E (hartree)
    double dE;           //!< |E_n - E_{n-1}|
    double commutator;   //!< [F,D] (the accelerator/DIIS error)
    double drho;         //!< relative charge-density change
    double order=0;      //!< the ORDER PARAMETER (tSCFIterator::SetOrderParameter); 0 when no probe is set
};

//! Everything one SCF step exposes to the per-system iteration display (doc/GPWPlan1.md item 2).  The base
//! \c Iterate builds a fully-populated trace each iteration -- all fields, all cheap -- and hands it to the
//! virtual \c DisplayColumns; each system's subclass renders only the columns it shows.  Real-valued
//! throughout (energies/residuals are real for both the rX and cX paths).  This bundle keeps the iterator's
//! members private: the display virtuals READ the trace, they do not reach back into the iterator.
struct IterationTrace
{
    size_t iteration;
    EnergyBreakdown eb;      //!< the full breakdown (Etotal, virial, per-term energies)
    double relax;            //!< the mixer step size α (the ρ_mix column's number)
    double relaxEff;         //!< α_eff: the step the mixer ACTUALLY delivered (== relax when unpreconditioned)
    const char* mixTag;      //!< ρ_mix identity: "Lin"/"Ker"/"Pul"
    const char* accelTag;    //!< accel identity: "Null"/"DIIS"/"GDM"
    int    accelCount;       //!< accel number (DIIS projection depth; 0 when not meaningful)
    double accelMinSV;       //!< DIIS conditioning: smallest singular value of the bordered B (NaN if none)
    double FD;               //!< [F,D] (the accelerator's orbital-gradient error)
    double dFD;              //!< Δ[F,D] = [F,D] − [F,D]_old
    double dRho;             //!< relative charge-density change ‖Δρ‖
    double dE;               //!< relative total-energy change ΔE/|E|
    size_t idealVirial;      //!< 2 (non-relativistic) or 1 (relativistic): the virial column shows GetVirial()+idealVirial
    bool   configChanged;    //!< did the occupied configuration change since the previous iteration? (the cfg '*')
    bool   lineSearch;       //!< this step ran the direct-min driver (GDM/OT): NO density mixing -> ρ_mix shows "----"
    double gridLostRel;      //!< signed grid-charge leak per electron, (∫ρ̃ − Tr(DS))/N (solids; 0 with no grid)
    // Frontier spectrum -- the gap column (solids show it always; molecules only under ReportBandGap()):
    double eHomo=0, eLumo=0, gap=0;
    bool   haveHomo=false, haveLumo=false, metallic=false, hole=false;
    // The caller's ORDER PARAMETER (tSCFIterator::SetOrderParameter), measured on THIS iteration's working
    // density.  nullptr name == no probe set == no column.
    const char* orderName=nullptr;   //!< the probe's short column label, e.g. "m_stag"
    double      order=0;             //!< its value this iteration
};

//! \brief ONE column of the SCF trace, carrying ALL THREE of its cells together.
//!
//! A trace column appears in three places -- the header label, the \c (tolerance) cell beneath it, and the
//! per-iteration value -- and they must agree.  They used to be written by three separate pieces of code,
//! so dropping a column meant remembering to drop its threshold too, and the lines silently stopped lining
//! up when you forgot.  Holding the three in one object makes that structural rather than a matter of care.
//!
//! A layout is therefore a LIST: \c AccumulateColumns builds it, and a subclass chooses the ORDER by when
//! it delegates to its base (V1.27, user design 2026-08-12).
struct ColumnData
{
    std::string title;                    //!< header label, e.g. "Δ[F,D]", "2+V/K", "ρ_lost/N"
    int         width = 10;               //!< DISPLAY width; padding is UTF-8 aware (Δ, ρ count as one column)
    //! The convergence threshold shown as "(1e-07)" under the label.  ABSENT (not a sentinel) for the
    //! columns that gate nothing -- ρ_lost/N is a health diagnostic, the gap and order parameter are probes.
    std::optional<double> tolerance;
    //! The per-iteration cell.  A callable rather than a field so the value comes from the SAME list as the
    //! label: a parallel "write the row" path is exactly what let the two drift apart before.
    std::function<void(std::ostream&, const IterationTrace&)> value;
};

// Templated on the matrix element type T (rX/cX); SCFIterator is the <double> alias (atoms/
// molecules), cSCFIterator the <dcmplx> instantiation that drives single-k plane-wave DFT through
// the same assemble->diagonalize->fill->build-density loop.
template <class T> class tSCFIterator
{
    typedef qchem::Hamiltonian::tHamiltonian<T>  ham_t;
    typedef qchem::WaveFunction::tWaveFunction<T>    wf_t;
    typedef qchem::WaveFunction::tSCFWaveFunction<T> scfwf_t;
    typedef qchem::SCFAccelerators::SCFAccelerator acc_t;   // NON-template manager (RealComplexPlan §6)
public:
    // The seed density is chosen by strategy (see ChargeDensity::SeedStrategy): Default resolves to
    // each path's present-day behaviour -- molecular -> CoreGuess, plane-wave -> Uniform.  \a st (the
    // structure) is needed only by the SAD seeds (atom Z + positions); null is fine otherwise.
    // \a basisOrtho / \a basisOrthoTol select how the orbital-overlap S is orthogonalised: Cholesky
    // (default; needs S positive-definite) or Eigen/SVD with a cutoff that drops near-null eigen/singular
    // values -- canonical orthogonalisation for a linearly-dependent basis (e.g. GPW's diffuse Gaussians
    // on a dense lattice, whose Bloch overlap goes singular).  \a basisOrthoTol<=0 keeps all.
    tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*,acc_t*,
                 ChargeDensity::SeedStrategy seed=ChargeDensity::SeedStrategy::Default,
                 const Structure* st=nullptr,
                 qchem::Ortho basisOrtho=qchem::Auto, double basisOrthoTol=0.0);
    //! Grid-continuation / explicit-seed ctor (doc/GPWPlan §0e): seed the SCF from a PRE-BUILT density
    //! \a seedDensity (TAKES OWNERSHIP -- consumed building the iteration-0 Fock, then freed) instead of
    //! resolving a SeedStrategy enum.  The intended use is grid continuation: a density converged on a
    //! CHEAP coarse collocation grid seeds an expensive fine-grid run on the SAME orbital basis, so the
    //! fine SCF starts in the physical basin rather than wandering into the −39 density/grid basin.  The
    //! orbital one-body matrices are grid-independent (only the collocation grid differs), so the coarse
    //! density matrix transfers directly, with no re-projection.  \a seedDensity's own basis block must
    //! outlive this ctor (it is read once here, in Init).
    tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*,acc_t*,
                 tChargeDensity<T>* seedDensity,
                 const Structure* st=nullptr,
                 qchem::Ortho basisOrtho=qchem::Auto, double basisOrthoTol=0.0);
    virtual ~tSCFIterator();
    virtual bool Iterate(const SCFParams& ipar);

    // Watch convergence live: the observer (if set) fires once per SCF iteration with
    // the current SCFProgress.  Read-only telemetry -- the observer must not drive the loop.
    using Observer = std::function<void(const SCFProgress&)>;
    void SetObserver(Observer obs) {itsObserver=std::move(obs);}

    //! Watch an ORDER PARAMETER die (doc/SymmetryUpgradePlan.md §9 "diagnostic metric"): a caller-supplied
    //! named scalar measured on the WORKING density every iteration, shown as a trace column and carried in
    //! SCFProgress.  A symmetry-broken solution (an AFM staggered moment, a charge disproportionation) is a
    //! basin the SCF can silently fall out of -- the seed is provably ordered, the answer provably isn't, and
    //! the converged numbers say nothing about WHICH iteration lost it.  The probe is the instrument that
    //! brackets it.  Deliberately a caller-supplied functor: the order parameter is system knowledge (which
    //! sites, which sign pattern), not something the iterator can infer.  Read-only telemetry, like the
    //! observer -- it must not touch the density.  \a name is a short column label ("m_stag"); an empty name
    //! or a null probe disables the column (the default: zero cost).
    using OrderProbe = std::function<double(const tDM_CD<T>&)>;
    void SetOrderParameter(const std::string& name, OrderProbe probe)
    {itsOrderName=name; itsOrderProbe=std::move(probe);}

    // SCFIterator drives the mutable SCFWaveFunction, but only ever hands clients the const
    // read view (they can query the converged state, never drive someone else's SCF loop).
    const wf_t* GetWaveFunction() const {return itsWaveFunction;}
    //! Grid-continuation MOM (doc/GPWPlan §0e): adopt \a from's converged occupied subspace as this run's FIXED
    //! MOM reference.  Call AFTER construction and BEFORE Iterate; takes effect with SCFParams::UseMOM (the
    //! reference is then held from iteration 1, never re-captured).  \a from must be a converged WF on the SAME
    //! orbital basis -- e.g. a coarse-density-grid solution seeding an expensive fine-grid run of the same
    //! system.  Seeding the DENSITY alone is not enough: on the fine grid a giant-response diffuse virtual sits
    //! at the frontier even at the physical density, so the occupied-subspace reference must transfer too.
    void AdoptMOMReference(const wf_t& from)
    {
        for (const auto& q : itsWaveFunction->GetQNs())
        {
            // Cross-cast to the policy's OrbitalView (abstract->abstract; TOrbitals implements it).
            auto* v=dynamic_cast<const qchem::OrbitalView<T>*>(from.GetOrbitals(q));
            assert(v && "AdoptMOMReference: the source orbitals must implement OrbitalView");
            itsOccState.AdoptReference(q, *v);   // R2.21: a reference is the STATE's memory, not the policy's
        }
    }
    EnergyBreakdown     GetEnergy() const;
    size_t              GetIterationCount() const {return itsIterationCount;}
    bool                Converged() const {return itsConverged;}

protected:
    // --- PER-SYSTEM iteration display (doc/GPWPlan1.md item 2) -------------------------------------------
    //! The verbose per-iteration trace, virtual so each system type presents the honest column set.  The
    //! base default is the MOLECULAR layout ({#, Etotal, [F,D], Δ[F,D], Δρ, 2+V/K, ρ_mix, accel, cfg}, plus
    //! the optional ReportBandGap() gap block); \c SolidSCFIterator overrides both for the solid/PP layout
    //! ({#, Etotal, [F,D], ΔE, Δρ, ρ_mix, accel, cfg, gap} -- no virial, gap always on).  \c Iterate calls
    //! these; the subclass never touches iterator internals (it reads only the \c IterationTrace / SCFParams).
    //! Build this layout's VARIABLE columns, in display order.  The base contributes the convergence gate
    //! (Δ[F,D]), Δρ, and -- iff the Hamiltonian says the virial theorem is meaningful (V1.27
    //! IsVirialValid) -- the virial.  A subclass calls the base and then swaps/inserts: \c Solid_ replaces
    //! the gate with ΔE/E and adds the collocation ρ_lost/N.  The INVARIANT frame around them (#, Etotal,
    //! [F,D], ρ_mix/accel/svMin/cfg, gap, order) is still written by the fixed writers below -- those cells
    //! carry their own alignment and flag glyphs, so they are not list-shaped yet.
    virtual void AccumulateColumns(std::vector<ColumnData>&, const SCFParams&, size_t idealVirial) const;
    virtual void DisplayColumnHeaders(std::ostream&, const SCFParams&, size_t idealVirial) const;
    virtual void DisplayColumns      (std::ostream&, const SCFParams&, const IterationTrace&) const;
    //! Is the frontier gap a PERMANENT column, or an optional diagnostic behind ReportBandGap()?
    //! PERIODIC-driven, not PP-driven: near-gapless flapping is a solid pathology.  (It travels with the
    //! grid columns today only because every gridded run here is also periodic -- see AccumulateColumns.)
    virtual bool GapIsPermanent() const {return false;}
    // Shared column writers the molecular/solid layouts compose (kept here so both subclasses reuse them;
    // the header-cell writers are file-local free functions in the .C, since they need no iterator state):
    void WriteRowPrefix   (std::ostream&, const IterationTrace&) const; //!< row:    #, Etotal, [F,D]
    void WriteMixAccelCfg (std::ostream&, const IterationTrace&) const; //!< row:    ρ_mix, accel, cfg
    void WriteGapColumn   (std::ostream&, const IterationTrace&) const; //!< row:    the frontier gap (+flags)
    void WriteOrderColumn (std::ostream&, const IterationTrace&) const; //!< row:    the order parameter (if probed)
    //! Header cell for the order-parameter column -- the label the caller gave SetOrderParameter (nothing
    //! when no probe is set), so both layouts announce the extra column the same way.
    void WriteHeadOrder   (std::ostream&) const;
    //! The Hamiltonian this run assembles, for the DISPLAY only: a layout has to ask it whether the virial
    //! theorem is meaningful (V1.27 IsVirialValid), and a subclass layout lives outside the base's privates.
    //! Const reference, not the pointer -- a layout may ASK the Hamiltonian, never re-seat it.
    const ham_t& Hamiltonian_() const {return *itsHamiltonian;}

    // --- PER-SYSTEM density mixing (doc/CleanupCandidates.md V1.10b) ------------------------------------
    //! \brief Build this run's density mixer.  Called once at the top of \c Iterate, after the seed density
    //! exists.  The base answer is the structure-neutral LINEAR D-mixer; \c SolidSCFIterator overrides with
    //! the periodic Kerker/Pulay G-space mixer when \c SCFParams asks for it.
    //!
    //! This is a virtual rather than a runtime probe on the geometry: "periodic or molecular?" is known by
    //! the class, so asking it again at run time (the old \c MakeDensityMixer capability probe, which fell
    //! back to linear mixing with a warning) put the decision one layer too low.  Everything the override
    //! needs is passed in, so iterator state stays private.
    virtual std::unique_ptr<qchem::ChargeDensity::tDensityMixer<T>>
        CreateMixer(const SCFParams& ipar, const tbs_t<T>* bs, const Structure* cell,
                    const tDM_CD<T>* seed) const;

private:
    typedef std::shared_ptr<tDM_CD<T>> cd_t;   //!< std-managed WORKING density (matrix-backed); no manual delete
    //! Seed the SCF: build the iteration-0 Fock from \a seed (a DFT-face tChargeDensity -- may be a fit, e.g.
    //! the SAD seed; or null for the core guess), diagonalize, and take the first real (matrix) density.
    //! \a bs / \a st are forwarded for the HF/DHF bootstrap (build a DFT sibling when the seed has no matrix
    //! but the Hamiltonian needs one -- see project_numericcd_refactor); null is fine for the core guess.
    void Initialize(tChargeDensity<T>* seed, const tbs_t<T>* bs, const Structure* st);
    //! The Hamiltonian's total energy for \a cd, with the wavefunction's Mermin −TS stamped into
    //! EnergyBreakdown::MinusTS so GetTotalEnergy() reads the free energy A=E−TS under Fermi smearing
    //! (0 with no smearing => plain E).  The entropy is a WF-side scalar (the Hamiltonian terms never
    //! see occupations), so the iterator -- which owns both -- is the seam that joins them.  doc/GPWPlan1.md 4b.
    EnergyBreakdown TotalEnergy(const tDM_CD<T>* cd) const;
    cd_t DirectMinStep(double Ecur, double mergeTol); //one direct-min step (returns new density; used by DirectMinDriver)

    void DisplayEigen   () const;

    //Raw ptrs owned, see destructor; the charge densities are std-managed (cd_t).
    ham_t*          itsHamiltonian;
    acc_t*          itsAccelerator;
    scfwf_t*        itsWaveFunction;
    //! The occupation SLOT (V1.11 inc 3; SCFStrategyPlan §6), now a STATE + a POLICY over it (R2.21).
    //!
    //! The STATE is the run's persistent memory -- MOM references, fill clocks, the −TS aggregate -- and is
    //! declared FIRST so it outlives every policy built over it.  It survives construction,
    //! AdoptMOMReference, reconfiguration and staged Iterate calls; and because it is scalar-INDEPENDENT,
    //! a mixed real/complex mesh keeps ONE ledger across blocks of both scalars (the mixed-mesh MOM case).
    //! The POLICY is the decision function, built here in its default (seed) state -- prescribed integer
    //! fill, kT=0, no MOM -- and configured from SCFParams at the top of Iterate, like the density mixer.
    qchem::OccupationState                      itsOccState;
    std::unique_ptr<qchem::OccupationPolicy<T>> itsOccPolicy
        = qchem::MakeOccupationPolicy<T>(qchem::OccupationConfig{}, itsOccState);
    cd_t            itsCD;       //!< current charge density (shared_ptr: lifetime by std, no reuse)
    cd_t            itsOldCD;    //!< previous charge density
    //! The SCF density lineage (ChargeDensity::Lineage).  SetWorkingCD makes each new itsCD the head, so a
    //! superseded density (itsOldCD, a stale copy) reports isActive()==false and trips the Hamiltonian's
    //! assert if reused.  Created in Init; one per SCF run.
    qchem::ChargeDensity::LineagePtr itsLineage;
    //! Assign the new working density AND make it the head of itsLineage (so the previous one goes inactive).
    //! Every itsCD (re)assignment in the working loop goes through here.
    void SetWorkingCD(cd_t cd) {itsCD=std::move(cd); if (itsCD && itsLineage) itsCD->JoinLineage(itsLineage);}

    size_t          itsIterationCount;
    bool            itsConverged;
    Observer        itsObserver;   //!< optional live-progress sink (default empty)
    OrderProbe      itsOrderProbe; //!< optional order-parameter probe (default empty == no column, no cost)
    std::string     itsOrderName;  //!< its column label

    // Density-face state.  The mixer (Linear / Kerker; see qchem.ChargeDensity.DensityMixer) is built per-run
    // from SCFParams at the top of Iterate and owns the mixing policy + state (relax, the Kerker ρ̃, ...).
    const tbs_t<T>*  itsBS = nullptr;         //!< orbital basis (for the Kerker G-space fit basis) -- from the ctor
    //! A PERSISTENT copy of the periodic cell (reciprocal lattice + volume for Kerker).  The ctor's raw \c st
    //! comes from a temporary (\c Lattice_3D::GetStructure returns a fresh \c make_shared), so it dangles by the
    //! time \c Iterate runs -- we deep-copy it here (periodic path only) so the Kerker mixer has a live cell.
    std::shared_ptr<const Structure> itsKerkerCell;
    std::unique_ptr<qchem::ChargeDensity::tDensityMixer<T>> itsMixer;  //!< the density-face concrete for this run
    // The two loop-face concretes (stateless).  Iterate selects one per macro-iteration by the accelerator's
    // WantsLineSearch() and dispatches Step() -- virtual dispatch in place of the old mode `if`.
    FixedPointDriver<T> itsFixedDriver;
    DirectMinDriver<T>  itsDirectDriver;
};

using SCFIterator  = tSCFIterator<double>;
using cSCFIterator = tSCFIterator<dcmplx>;

//! The ATOM/MOLECULE SCF driver (doc/GPWPlan1.md item 2).  Its display is the base default -- the virial-
//! bearing molecular column set -- so this subclass exists to NAME the molecular path (the Calculation /
//! AtomCalculation facades build it) and to give any future molecule-only columns a home.  Inherits every
//! base constructor unchanged.
class MolecularSCFIterator : public tSCFIterator<double>
{
public:
    using tSCFIterator<double>::tSCFIterator;
};

//! The SOLID / pseudopotential (GPW / plane-wave) SCF driver (doc/GPWPlan1.md item 2).  Overrides the trace
//! to the solid column set: the total-energy change ΔE gates it (a non-variational collocation SCF limit-
//! cycles [F,D]/Δρ above the fit floor while E settles), the virial is DROPPED (GTH local + KB projectors
//! break the Coulombic-homogeneity assumption behind 2+V/K), and the frontier gap is a PERMANENT column
//! (absorbing the former process-wide ReportBandGap() instrument -- the near-gapless flapping it watches is
//! a solid pathology).  Inherits every base (single-k dcmplx) constructor unchanged.
class SolidSCFIterator : public tSCFIterator<dcmplx>
{
public:
    using tSCFIterator<dcmplx>::tSCFIterator;
protected:
    //! The solid column list: the base's, with the gate swapped to ΔE/E and ρ_lost/N inserted.
    void AccumulateColumns(std::vector<ColumnData>&, const SCFParams&, size_t idealVirial) const override;
    bool GapIsPermanent() const override {return true;}   //!< near-gapless flapping is a solid pathology
    //! The periodic Kerker/Pulay G-space mixer when SCFParams asks for it (KerkerG0>0 or PulayDepth>0),
    //! else the base linear D-mixer.  A solid run HAS the periodic basis/cell/density by construction, so
    //! MakePeriodicMixer treats them as preconditions -- no capability probe, no silent fallback.
    std::unique_ptr<qchem::ChargeDensity::tDensityMixer<dcmplx>>
        CreateMixer(const SCFParams&, const tbs_t<dcmplx>*, const Structure*,
                    const tDM_CD<dcmplx>*) const override;
};

} //namespace


