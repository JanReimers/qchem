// File: SCFIterator/SCFIterator.C  Interface for an object that manages SCF convergence.
module;
#include <cassert>   // AdoptMOMReference's OrbitalView cross-cast guard
#include <memory>
#include <functional>
#include <string>
#include <iosfwd>
#include <vector>
#include <optional>
#include <string_view>
export module qchem.SCFIterator;
import qchem.SCFIterator.Types;
export import qchem.SCFAccelerator;
export import qchem.WaveFunction;
import qchem.WaveFunction.SCF;
import qchem.ElectronConfiguration.OccupationPolicy;   // the iterator's occupation SLOT (V1.11 inc 3)
export import qchem.SCFParams;
export import qchem.ChargeDensity.Seed;   // SeedStrategy / MakeSeedDensity
import qchem.LASolver;   // qchem::Ortho (the basis-overlap orthogonalisation knob, forwarded to the WF)
import qchem.BasisSet.Orbital_DFT_IBS;   // BasisSet::FIT_SF_ABS<T> (the G-space fit basis for Kerker rho-tilde extraction)
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

//! R2.7 (doc/SCFRobustnessPlan.md P2): NO DEAD GATE.  At the iteration that declares convergence, a LIVE [F,D] gate (MinFD set,
//! not off) reading EXACTLY 0.0 was never computed -- [F,D] arrives from the accelerator's GetError(), and a commutator
//! FDS-SDF is never exactly zero in floating point -- so the pass is the gate's absence, not its verdict (evidence G3: the Null
//! accelerator's [F,D] = 0 passed MinFD on every Kerker/Pulay run for months).  Throws.  Exported for the unit test.
//! NOT applied to Δρ/N: an exact 0 there is legitimate -- a density-INDEPENDENT Hamiltonian (a one-electron atom, A_*_E1/DE1)
//! reproduces its seed bit for bit at iteration 1 (measured 2026-10-30: 13 atom tests).
void RequireLiveGates(double FD, const Gate& minFD, size_t iteration);

//! \brief ONE CELL of the SCF iteration row: one column of the trace, as the console prints it AND as a plot wants it.
struct IterationCell
{
    enum class Scale { Linear, Log };      //!< how a plot shows it: Log for a convergence measure (its |value|)
    std::string           title;           //!< the header label ("Etotal", "[F,D]", "Δρ", "ρ_mix", "m_site", ...)
    std::optional<double> tolerance;       //!< the gate shown as "(1e-07)" under the label; absent for a column that gates nothing
    double                value;           //!< the number plotted -- NaN when the cell has none (a "----", a blank cfg)
    std::string           text;            //!< the cell exactly as the console prints it, padding stripped (flags, tags kept)
    Scale                 scale = Scale::Linear;
};

//! \brief ONE SCF ITERATION AS THE TRACE SHOWS IT (2026-10-31, OpenWork "The SCF iteration row has two homes"): the cells of
//! THIS run's layout, in display order -- the console writes them, the Observer receives them, built from ONE column list,
//! so a column added to a layout appears in both with no second edit.  It retired \c SCFProgress, whose five numbers were
//! a hand-picked subset (no Δ[F,D], ρ_lost/N, ρ_mix, accel, gap, ... and no titles or tolerances) with an ABSOLUTE |ΔE| that
//! the solid console printed RELATIVE under the same name; the row is now exactly the console's (controller ruling).
//! Built every iteration of a run that prints a trace OR has an observer.
struct IterationRow
{
    size_t                     iteration;
    std::vector<IterationCell> cells;      //!< excludes "#": this layout's columns, in display order
    //! \brief THIS iteration's full energy breakdown, of which the "Etotal" cell is \c GetTotalEnergy().
    //!
    //! Carried because a total energy alone cannot tell a healthy run from a collapsed one, but its TERMS
    //! can: a low-G charge runaway roughly DOUBLES the Hartree term \c Eee while the total stays a
    //! plausible negative number (MnO 2026-08-25: 13.48 -> 29.0 -> 35.1 Ha across three collapse arms,
    //! ranked correctly by \c Eee alone) -- and the probes read terms (Eee, the site moments) that are not columns.
    EnergyBreakdown            eb;
    //! The value of the cell titled \a title, NaN when this layout has no such column -- for a caller that consumes one
    //! number by NAME (the facade's detectors read "[F,D]" and "Δρ").
    double ValueOf(std::string_view title) const;
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
    const char* mixTag;      //!< ρ_mix identity: "Lin"/"Ker"/"Pul"
    const char* accelTag;    //!< the accelerator that PRODUCED this iterate ("DIIS"/"GDM"), "----" when none did (a plain or mixed step)
    int    accelCount;       //!< accel number (DIIS projection depth; 0 when not meaningful)
    double accelMinSV;       //!< DIIS conditioning: smallest singular value of the bordered B (NaN if none)
    double FD;               //!< [F,D] (the accelerator's orbital-gradient error)
    double dFD;              //!< Δ[F,D] = [F,D] − [F,D]_old
    double dRho;             //!< relative charge-density change ‖Δρ‖
    double dE;               //!< relative total-energy change ΔE/|E|
    size_t idealVirial;      //!< 2 (non-relativistic) or 1 (relativistic): the virial column shows GetVirial()+idealVirial
    bool   configChanged;    //!< did the occupied configuration change since the previous iteration? (the cfg '*')
    bool   lineSearch;       //!< a committed GEODESIC step (GDM/OT) produced this iterate: NO density mixing -> ρ_mix shows "----"
                             //!< (a direct-min iteration that fell back to a mixed step is false: it mixed)
    double gridLostRel;      //!< signed grid-charge leak per electron, (∫ρ̃ − Tr(DS))/N (solids; 0 with no grid)
    bool   gFilter;          //!< the mixer's step applies a G-dependent filter (tDensityMixer::HasGFilter; the XC ρ column)
    // Frontier spectrum -- the gap column (solids show it always; molecules only under ReportBandGap()):
    double eHomo=0, eLumo=0, gap=0;
    bool   haveHomo=false, haveLumo=false, metallic=false, hole=false;
    // The ORDER PARAMETER column (the signed max-|μ| integrated site moment): shown on every polarized run; "----" when the
    // run has no basins to integrate over.
    bool   hasOrder=false;           //!< this iteration's breakdown carried site moments
    double order=0;                  //!< the signed max-|μ| site moment
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
    //! The cell's NUMBER (the \c IterationCell::value an observer plots), from the same trace -- NaN-returning where the text
    //! is a "----".  Beside the text writer, so the number cannot come from a parallel path either.
    std::function<double(const IterationTrace&)> number;
    IterationCell::Scale scale = IterationCell::Scale::Linear;   //!< Log for a convergence measure
};

// Templated on the matrix element type T (rX/cX); SCFIterator is the <double> alias (atoms/
// molecules), cSCFIterator the <dcmplx> instantiation that drives single-k plane-wave DFT through
// the same assemble->diagonalize->fill->build-density loop.
//! What an explicit seed IS to the iterator (OpenWork §3 "A carried density is re-filled at kT = 0", 2026-10-21).  The
//! caller says it: only the composition root knows whether a density is a converged state of THIS graph -- a matrix-backed
//! seed is not enough (the Uniform and plane-wave strategy seeds are matrix-backed, in other block layouts).
enum class SeedRole
{
    Guess,           //!< build the iteration-0 Fock from it; the first working density is that Fock's fill (every strategy seed)
    WorkingDensity,  //!< a CONVERGED density of this basis and block layout (a restored state, an anneal / +-alpha / U-loop
                     //!< stage): when the run's occupation policy is not the default (kT > 0 or MOM) its iteration-0 fill is
                     //!< REDONE under that policy at the top of the first Iterate -- so a restored fixed point refills to itself
    ExactResume,     //!< a restored state with NO difference from this run (no WARM fingerprint entry): as WorkingDensity, and the
                     //!< carried D then IS the first working density -- the first Fock is built from the saved state itself, not
                     //!< from a diagonalization of it.  (OpenWork §3, 2026-10-29: as WorkingDensity, beta-MnO2's exact resume took
                     //!< two undamped steps before the mixer saw it -- the iteration-0 fill and iteration 1 -- and a 6e-6 residual
                     //!< came back as drho 1e-2; the fresh mixer then walked 0.25 Ha away.)  Not for a PERTURBED stage (+-alpha, a
                     //!< new U, an anneal): those need that first step (see Initialize).
};

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
                 qchem::Ortho basisOrtho=qchem::Auto, double basisOrthoTol=0.0, size_t numBands=0);
    //! Grid-continuation / explicit-seed ctor (doc/GPWPlan §0e): seed the SCF from a PRE-BUILT density
    //! \a seedDensity (TAKES OWNERSHIP -- consumed building the iteration-0 Fock, then freed) instead of
    //! resolving a SeedStrategy enum.  The intended use is grid continuation: a density converged on a
    //! CHEAP coarse collocation grid seeds an expensive fine-grid run on the SAME orbital basis, so the
    //! fine SCF starts in the physical basin rather than wandering into the −39 density/grid basin.  The
    //! orbital one-body matrices are grid-independent (only the collocation grid differs), so the coarse
    //! density matrix transfers directly, with no re-projection.  \a seedDensity's own basis block must
    //! outlive this ctor (it is read once here, in Init).
    //! \a role: \c SeedRole::WorkingDensity when \a seedDensity is a converged state of this graph (THROWS if it is not
    //! matrix-backed); the default \c Guess is the old behaviour.
    tSCFIterator(const tbs_t<T>*, const ElectronConfiguration*, ham_t*,acc_t*,
                 tChargeDensity<T>* seedDensity,
                 const Structure* st=nullptr,
                 qchem::Ortho basisOrtho=qchem::Auto, double basisOrthoTol=0.0, SeedRole role=SeedRole::Guess,
                 size_t numBands=0);   //!< \a numBands: lowest eigenpairs kept per block, 0 = all (LASolver::Factory)
    virtual ~tSCFIterator();
    //! Iterate resolves \c ipar.Mixer first (\c ResolveMixing), so \c Mixing::Auto becomes this run's mixer.
    virtual bool Iterate(const SCFParams& ipar);
    //! \a p with \c Mixer resolved for THIS run (never \c Auto): the mixer \c Iterate will build.  Public for the
    //! facade, which prints the mixer on its banner and caps a restart's LINEAR step -- both must see the resolved
    //! mixer, not the caller's \c Auto.
    SCFParams ResolveMixing(const SCFParams& p) const;
    //! R1.3: print \a h in the one hand-off format; \a reseat re-seats the field mixer on the current D (history dropped).
    //! Returns whether anything was re-seated.
    bool HandOff(const HandOffRequest& h, bool reseat);

    // Watch convergence live: the observer (if set) fires once per SCF iteration with that iteration's
    // IterationRow -- the console's own cells.  Read-only telemetry -- the observer must not drive the loop.
    using Observer = std::function<void(const IterationRow&)>;
    void SetObserver(Observer obs) {itsObserver=std::move(obs);}

    // (No SetOrderParameter any more -- 2026-09-14.  The order parameter is the INTEGRATED site moment the
    //  spin-native XC term writes into every iteration's EnergyBreakdown (R1.0h); the m_site column reads it
    //  from there.  A caller-supplied point probe cannot be a moment.)

    // SCFIterator drives the mutable SCFWaveFunction, but only ever hands clients the const
    // read view (they can query the converged state, never drive someone else's SCF loop).
    const wf_t* GetWaveFunction() const {return itsWaveFunction.get();}
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
    //! \brief Why the last \c Iterate STOPPED on a physically impossible first energy instead of iterating (empty = it did not).
    //! The first iterate, the seed's Hamiltonian diagonalized once, lies close to the answer on any working Hamiltonian, so a first
    //! energy past the layout's bound (\c FirstEnergyRefusal) names a broken Hamiltonian (2026-11-01, beta-MnO2 on a dense-d basis:
    //! +4.19e4 Ha, an aliased raster amplified by an ill-conditioned S), and iterating it only turns the defect into a
    //! plausible-looking non-convergence.  The facade turns it into its failed Outcome.
    const std::string&  UnphysicalStop() const {return itsUnphysicalStop;}

protected:
    // --- PER-SYSTEM iteration display (doc/GPWPlan1.md item 2) -------------------------------------------
    //! The verbose per-iteration trace, virtual so each system type presents the honest column set.  The
    //! base default is the MOLECULAR layout ({#, Etotal, [F,D], ΔE/E, [Δ[F,D]], Δρ, 2+V/K, ρ_mix, accel, cfg}, plus
    //! the optional ReportBandGap() gap block); \c SolidSCFIterator overrides both for the solid/PP layout
    //! ({#, Etotal, [F,D], ΔE/E, Δρ, ρ_lost/N, XC ρ, ρ_mix, accel, cfg, gap} -- gap always on).  \c Iterate calls
    //! these; the subclass never touches iterator internals (it reads only the \c IterationTrace / SCFParams).
    //! Build this layout's VARIABLE columns, in display order.  The base contributes ΔE/E, Δ[F,D] only when MinΔFD gates,
    //! Δρ, and -- iff the Hamiltonian says the virial theorem is meaningful (V1.27
    //! IsVirialValid) -- the virial.  A subclass calls the base and then inserts: \c Solid_ adds the collocation ρ_lost/N
    //! and XC ρ after Δρ.  The INVARIANT frame around them (Etotal, [F,D],
    //! ρ_mix / accel / svMin / cfg, gap, m_site) joins them in \c RowColumns (2026-10-31), each cell carrying its own
    //! alignment and flag glyphs inside its text writer.
    virtual void AccumulateColumns(std::vector<ColumnData>&, const SCFParams&, size_t idealVirial) const;
    //! THE WHOLE ROW, in display order ("#" aside): the invariant frame around this layout's \c AccumulateColumns.  The ONE
    //! list both the console (\c DisplayColumns) and the observer (\c MakeRow) are built from.
    std::vector<ColumnData> RowColumns(const SCFParams&, size_t idealVirial) const;
    //! This iteration's row for the observer: every \c RowColumns cell's title, tolerance, number, text and scale.
    IterationRow MakeRow(const SCFParams&, const IterationTrace&) const;
    virtual void DisplayColumnHeaders(std::ostream&, const SCFParams&, size_t idealVirial) const;
    virtual void DisplayColumns      (std::ostream&, const SCFParams&, const IterationTrace&) const;
    //! Is the frontier gap a PERMANENT column, or an optional diagnostic behind ReportBandGap()?
    //! PERIODIC-driven, not PP-driven: near-gapless flapping is a solid pathology.  (It travels with the
    //! grid columns today only because every gridded run here is also periodic -- see AccumulateColumns.)
    virtual bool GapIsPermanent() const {return false;}
    //! The iteration-1 refusal (\c UnphysicalStop): the reason \a eb (\a nElectrons in the system) is impossible for this layout's
    //! system, or empty.  A FINITE system's energy measured from the electrons dispersed to infinity, E minus the structure constants
    //! (rest mass and -TS included), is <= 0 for a ground state at any geometry and charge (the TOTAL can be positive: H2 at
    //! R = 0.1 bohr is +7 Ha) -- so a first iterate above 0 is refused.
    virtual std::string FirstEnergyRefusal(const EnergyBreakdown& eb, double nElectrons) const;
    // The frame cells are ColumnData now (RowColumns); the gap writer stays a member because the OPTIONAL frontier block
    // (ReportBandGap / QCHEM_DIAGNOSTICS=frontier) appends it after the row on a molecular layout, where it is not a column.
    void WriteGapColumn   (std::ostream&, const IterationTrace&) const; //!< row:    the frontier gap (+flags)
    //! Header cell for the order-parameter column ("m_site" on a polarized run, nothing otherwise), so
    //! both layouts announce the extra column the same way.
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
    //! Can this run mix the density in G space (Kerker / Pulay)?  Decides what \c ResolveMixing turns \c Mixing::Auto into.
    virtual bool MixesInGSpace(const tbs_t<T>*, const Structure*, const tDM_CD<T>*) const {return false;}

private:
    typedef std::shared_ptr<tDM_CD<T>> cd_t;   //!< std-managed WORKING density (matrix-backed); no manual delete
    //! Seed the SCF: build the iteration-0 Fock from \a seed (a DFT-face tChargeDensity -- may be a fit, e.g.
    //! the SAD seed; or null for the core guess), diagonalize, and take the first real (matrix) density.
    //! \a bs / \a st are forwarded for the HF/DHF bootstrap (build a DFT sibling when the seed has no matrix
    //! but the Hamiltonian needs one -- see project_numericcd_refactor); null is fine for the core guess.
    void Initialize(tChargeDensity<T>* seed, const tbs_t<T>* bs, const Structure* st, SeedRole role=SeedRole::Guess);
    //! The Hamiltonian's total energy for \a cd, with the wavefunction's Mermin −TS stamped into
    //! EnergyBreakdown::MinusTS so GetTotalEnergy() reads the free energy A=E−TS under Fermi smearing
    //! (0 with no smearing => plain E).  The entropy is a WF-side scalar (the Hamiltonian terms never
    //! see occupations), so the iterator -- which owns both -- is the seam that joins them.  doc/GPWPlan1.md 4b.
    EnergyBreakdown TotalEnergy(const tDM_CD<T>* cd) const;
    cd_t DirectMinStep(double Ecur, double mergeTol); //one direct-min step (returns new density; used by DirectMinDriver)

    void DisplayEigen   () const;

    //! \name The SCF's three collaborators, split by WHO MADE THEM (R2.22)
    //!
    //! The rule, and it is the whole of the ownership story here: **what this class is HANDED it does not
    //! delete; what this class MAKES it holds in a `unique_ptr`.**  Before R2.22 all three were raw and all
    //! three were deleted, which made the iterator destroy a Hamiltonian built by its own caller -- and
    //! THAT is why an anneal stage had to rebuild one (15.5 s/call, 19% of a threaded MnO run): the
    //! composition root could not keep the object it owned.
    //!@{
    //! NON-OWNING.  The composition root (SolidCalculation / MoleculeCalculation / AtomCalculation, or a test
    //! harness) owns the Hamiltonian and must outlive this iterator.  It is a pure function of
    //! (structure, basis, species, functional, xcMesh, vxcFit), so ONE serves a whole annealed run.
    ham_t*          itsHamiltonian;
    //! NON-OWNING, and deliberately REPLACED per anneal stage by the owner: a stage change invalidates the
    //! Pulay/DIIS history, and the accelerator TYPE itself changes (anneal on Ladder, finish on GDM).
    //! Cheap to construct, so rebuilding it is right -- it is the Hamiltonian that was never the problem.
    acc_t*          itsAccelerator;
    //! OWNED -- this one the iterator genuinely makes (\c WaveFunction::Factory, in the ctor's init list).
    std::unique_ptr<scfwf_t> itsWaveFunction;
    //!@}
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
    bool            itsGeodesicTaken=false;   //!< the last \c DirectMinStep committed a geodesic (else it fell back to a MIXED step)
    std::string     itsUnphysicalStop;   //!< \c UnphysicalStop: empty unless the last Iterate refused its first energy
    std::shared_ptr<tDM_CD<T>> itsCarried;   //!< a WorkingDensity / ExactResume seed, held to the first Iterate (see SeedRole)
    bool            itsCarriedIsWorking=false;   //!< ExactResume: the carried D becomes the working density there
    Observer        itsObserver;   //!< optional live-progress sink (default empty)

    // Density-face state.  The mixer (Linear / Kerker; see qchem.ChargeDensity.DensityMixer) is built per-run
    // from SCFParams at the top of Iterate and owns the mixing policy + state (relax, the Kerker ρ̃, ...).
    const tbs_t<T>*  itsBS = nullptr;         //!< orbital basis (for the Kerker G-space fit basis) -- from the ctor
    //! A PERSISTENT copy of the periodic cell (reciprocal lattice + volume for Kerker).  The ctor's raw \c st
    //! comes from a temporary (\c Lattice_3D::GetStructure returns a fresh \c make_shared), so it dangles by the
    //! time \c Iterate runs -- we deep-copy it here (periodic path only) so the Kerker mixer has a live cell.
    std::shared_ptr<const Structure> itsKerkerCell;
    //! D-MINDRHO: \f$\int|\rho_{out}-\rho_{in}|\,d^3r/N\f$ -- the Hamiltonian's own quadrature when a term owns one (the
    //! XC term), else a mesh built ONCE from a copy of the structure (HF, 1-electron).
    double ResidualPerElectron(const tChargeDensity<T>& in, const tChargeDensity<T>& out) const;
    std::shared_ptr<const Structure> itsResidualStructure;   //!< deep copy (the ctor's st may be a temporary)
    mutable rvec3vec_t itsResidualPts;                       //!< the fallback mesh, built on first need
    mutable rvec_t     itsResidualWts;
    double             itsResidual=0.0;                      //!< this iteration's measure (the driver writes it)
    std::unique_ptr<qchem::ChargeDensity::tDensityMixer<T>> itsMixer;  //!< the density-face concrete for this run
    qchem::ChargeDensity::tAdaptiveMixer<T>* itsAdaptive=nullptr;      //!< its step-adapting face, if it has one (cross-cast once)
    // The two loop-face concretes (stateless).  Iterate selects one per macro-iteration by the accelerator's
    // WantsLineSearch() and dispatches Step() -- virtual dispatch in place of the old mode `if`.
    FixedPointDriver<T> itsFixedDriver;
    DirectMinDriver<T>  itsDirectDriver;
};

using SCFIterator  = tSCFIterator<double>;
using cSCFIterator = tSCFIterator<dcmplx>;

//! The ATOM/MOLECULE SCF driver (doc/GPWPlan1.md item 2).  Its display is the base default -- the virial-
//! bearing molecular column set -- so this subclass exists to NAME the molecular path (the MoleculeCalculation /
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
    //! A periodic cell has NO electrons-at-infinity limit (Ewald: E_el +1.59 Ha on a healthy Si Γ GTH run), and even its TOTAL can be
    //! positive when healthy (Na q1 at the atomic density: +0.047 Ha), so there is no theorem here, only a stated SCALE: a first
    //! iterate above +1 Ha PER ELECTRON is refused.  Measured: healthy cells sit within a few Ha per electron below 0 to +0.05 above;
    //! the broken dense-d β-MnO₂ first iterates were +22 and +1100 Ha per electron (OpenWork row: a sounder periodic bound).
    std::string FirstEnergyRefusal(const EnergyBreakdown& eb, double nElectrons) const override;
    //! The resolved \c SCFParams::Mixer: the periodic Kerker or Kerker/Pulay G-space mixer, or the base linear D-mixer.
    //! A solid run HAS the periodic basis/cell/density by construction, so the G-space mixer factories treat them as
    //! preconditions -- no capability probe, no silent fallback.
    std::unique_ptr<qchem::ChargeDensity::tDensityMixer<dcmplx>>
        CreateMixer(const SCFParams&, const tbs_t<dcmplx>*, const Structure*,
                    const tDM_CD<dcmplx>*) const override;
    //! Yes when the G-space factories' precondition holds (basis, periodic cell, a density with a Fourier face) -- GPW and PW
    //! alike; an iterator built without its Structure (no cell) mixes linearly under Auto.
    bool MixesInGSpace(const tbs_t<dcmplx>*, const Structure*, const tDM_CD<dcmplx>*) const override;
};

} //namespace


