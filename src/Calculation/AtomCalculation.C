// File: Calculation/AtomCalculation.C
//
// qchem::AtomCalculation -- the production front door for a single-atom SCF (the atom/Dirac/PP analogue
// of qchem::MoleculeCalculation, which is molecule-only by deliberate design).  It owns the whole object graph
// (atom + atomic EC + atomic basis + Hamiltonian + accelerator + iterator), runs the assemble->converge
// recipe, and exposes high-level questions: Energy(), EnergyTerms(), Density(), Orbital(i).
//
// This is the single tested atom recipe that previously lived only in the QchemTester/TestAtom test
// scaffold.  A SIBLING of qchem::MoleculeCalculation rather than an extension: the molecular facade is being
// consumed by the GUI team and must stay stable, and atoms genuinely differ (atomic exponent-pool bases,
// atomic / Dirac / pseudo ECs, the degree-0 (1-direction) atomic mesh, relativistic + PP models the molecule facade
// has no business expressing).  The two share the SCFIterator converge machinery and AcceleratorOptions.
module;
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <optional>
#include <nlohmann/json.hpp>
export module qchem.AtomCalculation;

import qchem.Structure;             // Structure, Atom
import qchem.Types;                 // hmat_t (a restored block density)
import qchem.ScalarFunction;        // ScalarFunction<double>
import qchem.BasisSet;              // BasisSet::Real_BS
import qchem.BasisSet.Radial.Factory; // BasisSet::Radial::Type, BasisSetAccuracy (the atomic exponent-pool basis)
import qchem.Hamiltonian.Factory;   // Hamiltonian::Model, SpinGroup
import qchem.Mesh;                  // qcMesh::MeshParams
import qchem.ElectronConfiguration; // ElectronConfiguration
import qchem.SCFIterator;           // SCFIterator, SCFParams, IterationRow, EnergyBreakdown
import qchem.Symmetry.Irrep;        // Irrep
export import qchem.Calculation;    // the abstract base of the three facades
export import qchem.StateView;      // ViewState's answer: a saved state, drawn
export import qchem.Symmetry.Spin;  // Spin (per-spin irrep enumeration), SpinGroup (an RunOptions field)
import qchem.Orbitals;              // Orbital, Orbitals
import qchem.ChargeDensity;         // rDM_CD
import qchem.ChargeDensity.Seed;    // SeedStrategy
import qchem.MoleculeCalculation;           // reuse AcceleratorOptions (read-only; does NOT modify MoleculeCalculation)
export import qchem.RunOptions;             // RunOptions: the one options tree (V-FACADE-PARITY P0b-2)
export import qchem.AtomSetup;              // AtomSetup: the run setup an AtomCalculation is built from (V-FACADE-PARITY P1)

export namespace qchem
{

using Hamiltonian::Model;                 // {E1, HF, DE1, DHF, Xalpha, LDA}

class AtomCalculation : public virtual tCalculation<double>
{
public:
    using sf_t       = ScalarFunction<double>;
    using orbitals_t = qchem::Orbitals::Orbitals;
    using orbital_t  = qchem::Orbitals::Orbital;
    using Observer   = qchem::SCFIterator::SCFIterator::Observer;

    //! Build the whole graph from \a setup (\c AtomSetup::Resolve: Z, charge, options, accelerator) and converge ONCE with \a params.
    //! The EC is chosen from the model: Dirac models get an AtomDirac_EC, the rest an Atom_EC.  (Atom SCFs often want per-Z params and
    //! some are expensive, so the params go in here -- a single converge -- rather than a default-then-reconverge.)
    explicit AtomCalculation(std::shared_ptr<const AtomSetup> setup, const SCFParams& params = {});
    ~AtomCalculation();

    AtomCalculation(const AtomCalculation&)            = delete;
    AtomCalculation& operator=(const AtomCalculation&) = delete;

    //! \name BUILT, NOT CONVERGED (2026-11-02: the twin of \c MoleculeCalculation::Prepare / \c SolidCalculation::Prepare)
    //! The whole graph -- basis, Hamiltonian with every mesh, accelerator, SCF iterator (its wave-function build files
    //! basis.perIrrep) -- into the CALLER'S report run (a facade never opens one, user ruling D 2026-11-02); no SCF.  Run it with
    //! \c Converge(params), whose FIRST call reuses this graph.  Until then the accessors have no SCF to read.
    //!@{
    //! REFUSED (a value, R9) when a policy refuses while the graph is built (a 4-centre tensor that cannot fit in memory, a PP the run
    //! cannot hold); its policy row is already filed.
    static Outcome<std::unique_ptr<AtomCalculation>, BuildRefusal> Prepare(std::shared_ptr<const AtomSetup> setup, Observer onIteration = nullptr);
    //!@}
    //! \name RESTART FROM A SAVED STATE (V-FACADE-STATE; \c SaveState / \c JudgeState are the base's)
    //! The atomic twins of \c MoleculeCalculation::Restart / \c PrepareRestart, built from a setup as the constructor is:
    //! REFUSED on another kind, another radial basis, other electrons or blocks; an EXACT RESUME (the saved D drives the first Fock,
    //! a Ladder starts on its GDM rung at kT = 0) when nothing differs; else a WARM START with every difference named.
    //!@{
    static Outcome<std::unique_ptr<AtomCalculation>, RestartRefusal>
        Restart(const std::string& path, std::shared_ptr<const AtomSetup> setup, const SCFParams& params = {}, Observer onIteration = nullptr);
    //! \c Restart up to its SCF; run it with \c Converge.
    static Outcome<std::unique_ptr<AtomCalculation>, RestartRefusal>
        PrepareRestart(const std::string& path, std::shared_ptr<const AtomSetup> setup, const SCFParams& params = {},
                       Observer onIteration = nullptr);
    const std::optional<RestartSummary>& RestartedFrom() const override;   //!< \c PrepareRestart's judgement; empty for a fresh run
    //!@}

    //! \brief A SAVED STATE, DRAWN (\c StateView, lane G 2026-11-13): the state at \a path on the basis \a setup builds (as
    //! \c PrepareRestart builds it), then the saved D, C, eps, f -- no Hamiltonian, no SCF.  REFUSED as \c Calculation::JudgeView
    //! refuses, never thrown.  Reports one level in, under the scope its caller opened (a \c report::Section "view").
    static Outcome<std::unique_ptr<StateView>, RestartRefusal>
        ViewState(const std::string& path, std::shared_ptr<const AtomSetup> setup);

    //! \c Converge(edit) is the base's template method (V-FACADE-PARITY P3): the setup's \c opts.scf, edited by a one-field override
    //! (R10); the first call reuses a \c Prepare'd graph.  \c Result() (the base's) re-reads its outcome.
    void OnIteration(Observer obs);

    double                 Energy()      const;
    qchem::EnergyBreakdown EnergyTerms() const;
    double                 TotalCharge() const;

    const sf_t&       Density()         const;
    const sf_t&       HOMO()            const;
    const sf_t&       Orbital(size_t i) const;
    size_t            NumOccupied()     const;
    const orbitals_t* Orbitals(const Irrep&) const;

    //! Lower-level orbital access for the atomic/Dirac diagnostics (eigenvalues, kappa, fine structure):
    //! enumerate the per-spin irreps, and fetch the index-th orbital within one irrep.
    BasisSet::irrepv_t GetIrreps(const Spin& ms) const;
    const orbital_t*   GetOrbital(size_t index, const Irrep& qns) const;

    size_t           IterationCount() const;   //!< the run's SCF iterations: the entry SCF plus every atom-pilot re-convergence
    bool             IsConverged()    const;
    int              GetZ()           const {return itsSetup->Z;}
    const Structure& GetStructure()   const {return *itsSetup->structure;}
    const AtomSetup& GetSetup()       const override {return *itsSetup;}   //!< the setup this atom was built from (a restart or a sibling reads its decisions); the base's hook, covariant

private:
    typedef std::pair<double, const sf_t*> occ_t;
    struct BuildOnly {};
    //! Everything but the SCF: the members and the electron configuration (the graph is \c BuildGraph's).
    AtomCalculation(BuildOnly, std::shared_ptr<const AtomSetup> setup, Observer onIteration);
    //! The accelerator, the basis (first time), the Hamiltonian and the SCF iterator -- inside the caller's open report run.
    void BuildGraph();
    std::string AcceleratorTag() const;   //!< the accelerator type the setup resolved, as the scf header names it
    bool itsPrepared = false;             //!< \c Prepare built the graph: the next Converge reuses it instead of rebuilding
    void RebuildSampling();
    //! \name \c Converge's STEPS (\c tCalculation's hooks)
    //!@{
    //! One stage: the recipe's report, the graph (\c Prepare's, reused once, else rebuilt), the SCF, the grid pilot, the result line.
    bool RunStage(const SCFParams& params) override;
    //! A schedule's stage (P4): on the setup's accelerator and recipe (\c AdmitStageOnSetup refuses another), from the density the last
    //! stage reached once any SCF has run -- stage 0 of a \c Prepare'd graph runs on that graph, from the setup's seed.
    void EnterStage(const SCFStage& stage, size_t index) override;
    //!@}
    //! \name THE PROOF'S ANSWERS (\c tCalculation's hooks)
    //!@{
    std::optional<SCFFailure>             MakeFailure()      const override;
    EnergyBreakdown                       ProofEnergyTerms() const override;
    double                                ProofCharge()      const override;
    size_t                                ProofIterations()  const override;
    const ChargeDensity::tDM_CD<double>&  ProofDensity()     const override;
    const ScalarFunction<double>&         ProofSpinDensity() const override;
    //!@}
    mutable std::unique_ptr<const sf_t>   itsSpin;   //!< m(r) of a polarized attempt, built on first use (reset by \c RebuildSampling)
    //! \name THE STATE HOOKS (\c Calculation's template method, V-FACADE-STATE)
    //!@{
    void                    RecordState(StateRecord&, double kT) const override;   //!< the base's keys + the radial basis, model, mesh
    std::string             StateKind() const override {return "atom";}
    const Structure&        StateStructure() const override {return *itsSetup->structure;}
    std::string             StateFunctional() const override;
    std::string             KindTag() const override;
    double                  StateKT() const override;
    std::vector<StateBlock> StateBlocks(bool withOrbitals) const override;
    StateSummary            MakeStateSummary() const override;
    //!@}
    //! A RESTORED density, per block in the canonical order, consumed by the next \c BuildGraph (\c PrepareRestart).
    std::vector<hmat_t<double>>   itsRestored;
    std::unique_ptr<qchem::ChargeDensity::rDM_CD> RestoredDensity() const;   //!< \c itsRestored as a density on this basis
    std::string FunctionalName() const;                                        //!< the model as the fingerprint and summary name it
    bool                          itsRestoredExact = false;
    bool                          itsStartOnGDM    = false;   //!< the next accelerator is a Ladder started on its GDM rung, once
    double                        itsLastkT        = 0.0;
    std::optional<RestartSummary> itsRestartedFrom;

    std::shared_ptr<const AtomSetup>     itsSetup;       //!< every decision resolved once: Z, Ne, the options, the accelerator, the seed (its structure the Hamiltonian shares)
    ElectronConfiguration*               itsEC    = nullptr;  //!< owned
    BasisSet::Real_BS*                   itsBasis = nullptr;  //!< owned
    //! The run's Coulomb-fit slot (doc/IntegralCacheOwnershipPlan.md S2b): every Hamiltonian builder demands one.  A slot fills ONCE, so the
    //! atom pilot replaces it when it raises the Dunlap fit's rung (a unique_ptr for that reason, never for sharing).
    std::unique_ptr<Hamiltonian::CoulombFitSlot<double>> itsCoulombFit=std::make_unique<Hamiltonian::CoulombFitSlot<double>>();
    //! \brief The atom's RESOLVED grids (policy: \c ResolveGrids at construction, raised by \c PilotGrids): the XC radial count, the Hartree
    //! route with its radial count, and the Dunlap fit's rung.  "stated" = the caller gave it, so it is audited and never raised.
    struct Grids
    {
        int  xcRadial=50;
        bool xcStated=false;
        Hamiltonian::CoulombRoute coulomb;
        bool hStated=false;
        int  fitRung=0;
    };
    Grids itsGrids;
    bool  itsDFT=false;       //!< the model reads a density on a mesh (DFT or the PP pseudo-atom): only these have grids to size
    bool  itsPiloted=false;   //!< the atom pilot ran (once per calculation)
    size_t itsIterations=0;   //!< the last Converge's SCF iterations, the pilot's re-convergences included (the iterator counts only its own)
    std::unique_ptr<qchem::ChargeDensity::rDM_CD> itsCarried;   //!< the converged D the pilot carries into its re-convergence, consumed by \c BuildGraph
    void ResolveGrids();
    //! The electron configuration the SCF occupies: \c opts.configuration STATED (checked, refused when it cannot run) or derived; reported
    //! under policy AtomConfiguration when derived.
    ElectronConfiguration* MakeConfiguration() const;
    //! \return the SCF iterations its re-convergences on raised grids spent (0 when nothing was raised).
    size_t PilotGrids(const SCFParams& params);
    SCFParams PilotWarmStart(const SCFParams& params) const;   //!< the pilot re-convergence's SCF parameters (policy AtomPilotWarmStart)
    //! A Hamiltonian on \a g, its Coulomb fit from \a slot (filled here with \a g's rung when it is above 0).
    qchem::Hamiltonian::rHamiltonian* MakeHamiltonianOn(const Grids& g, const Hamiltonian::CoulombFitSlot<double>& slot) const;
    //! R2.22: the iterator no longer deletes what it is handed, so this facade owns the Hamiltonian and
    //! the accelerator.  Both are still REBUILT per Converge() -- a molecular Hamiltonian is cheap, so
    //! that stays as it was; assigning here frees the previous pair at exactly the old moment.  Declared
    //! BEFORE itsScf so they outlive the iterator built over them.
    std::unique_ptr<qchem::Hamiltonian::rHamiltonian>       itsHam;
    std::unique_ptr<qchem::SCFAccelerators::SCFAccelerator> itsAccel;
    qchem::SCFIterator::SCFIterator*     itsScf   = nullptr;  //!< owned (the Hamiltonian/accelerator are NOT its)
    std::unique_ptr<qchem::ChargeDensity::rDM_CD> itsDensity;
    std::vector<occ_t>                   itsOccupied;
    Observer                             itsObserver;
};

} // namespace qchem
