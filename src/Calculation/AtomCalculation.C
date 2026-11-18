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

export namespace qchem
{

using Hamiltonian::Model;                 // {E1, HF, DE1, DHF, Xalpha, LDA}

//! The atom's DEFAULT XC mesh (unstated \c RunOptions::xcMesh): the RADIAL-ONLY set (\c AngularKind::None, one direction at 4 pi), exact
//! because the atom's density is spherical by construction; \c nRadial 0 = unstated: the run enters at 50 MHL points and climbs to eps.
//! Was \c RunOptions::mesh's member default.
inline const qcMesh::MeshParams kAtomMeshDefault = {.radial = qcMesh::RadialKind::MHL, .nRadial = 0, .mhl_m = 3, .mhl_alpha = 2.0,
                                                    .angular = qcMesh::AngularKind::None, .angularDegree = 0, .beckeOrder = 2};
//! THE ATOM'S VIEW of the options tree: \a asked with the atom's defaults filled (model HF, \c kAtomMeshDefault, ortho Auto) and what an
//! atom does not consume REFUSED (policy ForeignOption): the molecule and solid blocks, Nelec / charge / multiplicity (the atom takes
//! Z and charge as arguments, its spin from the configuration), pp / species (the atom block's pseudopotential + valence), the
//! imposition, spins, the top-level accelerator / recipe, orthoTol, label, saveStateTo.
RunOptions AtomOptions(const RunOptions& asked);
//! A caller that edits ONE field of an atom's XC mesh starts from the atom's default (\c kAtomMeshDefault), filled on first use.
inline qcMesh::MeshParams& AtomXCMesh(RunOptions& o) { if (!o.xcMesh) o.xcMesh=kAtomMeshDefault; return *o.xcMesh; }

class AtomCalculation : public virtual Calculation
{
public:
    using sf_t       = ScalarFunction<double>;
    using orbitals_t = qchem::Orbitals::Orbitals;
    using orbital_t  = qchem::Orbitals::Orbital;
    using Observer   = qchem::SCFIterator::SCFIterator::Observer;

    //! Build the whole graph and converge ONCE with \a params.  \a Z is the nuclear charge, \a charge the
    //! net ionic charge (electrons = Z - charge).  The EC is chosen from the model: Dirac models get an
    //! AtomDirac_EC, the rest an Atom_EC.  (Atom SCFs often want per-Z params and some are expensive, so
    //! the params go in here -- a single converge -- rather than a default-then-reconverge.)
    explicit AtomCalculation(int Z, int charge = 0,
                             const RunOptions&    opts   = {},
                             const SCFParams&          params = {},
                             const AcceleratorOptions& acc    = {});
    ~AtomCalculation();

    AtomCalculation(const AtomCalculation&)            = delete;
    AtomCalculation& operator=(const AtomCalculation&) = delete;

    //! \name BUILT, NOT CONVERGED (2026-11-02: the twin of \c MoleculeCalculation::Prepare / \c SolidCalculation::Prepare)
    //! The whole graph -- basis, Hamiltonian with every mesh, accelerator, SCF iterator (its wave-function build files
    //! basis.perIrrep) -- into the CALLER'S report run (a facade never opens one, user ruling D 2026-11-02); no SCF.  Run it with
    //! \c Converge(params), whose FIRST call reuses this graph.  Until then the accessors have no SCF to read.
    //!@{
    static std::unique_ptr<AtomCalculation> Prepare(int Z, int charge = 0, const RunOptions& opts = {},
                                                    const AcceleratorOptions& acc = {}, Observer onIteration = nullptr);
    //!@}
    //! \name RESTART FROM A SAVED STATE (V-FACADE-STATE; \c SaveState / \c JudgeState are the base's)
    //! The atomic twins of \c MoleculeCalculation::Restart / \c PrepareRestart, built from the same arguments as the constructor:
    //! REFUSED on another kind, another radial basis, other electrons or blocks; an EXACT RESUME (the saved D drives the first Fock,
    //! a Ladder starts on its GDM rung at kT = 0) when nothing differs; else a WARM START with every difference named.
    //!@{
    static Outcome<std::unique_ptr<AtomCalculation>, RestartRefusal>
        Restart(const std::string& path, int Z, int charge = 0, const RunOptions& opts = {}, const SCFParams& params = {},
                const AcceleratorOptions& acc = {}, Observer onIteration = nullptr);
    //! \c Restart up to its SCF; run it with \c Converge.
    static Outcome<std::unique_ptr<AtomCalculation>, RestartRefusal>
        PrepareRestart(const std::string& path, int Z, int charge = 0, const RunOptions& opts = {}, const SCFParams& params = {},
                       const AcceleratorOptions& acc = {}, Observer onIteration = nullptr);
    const std::optional<RestartSummary>& RestartedFrom() const override;   //!< \c PrepareRestart's judgement; empty for a fresh run
    //!@}

    //! \brief A SAVED STATE, DRAWN (\c StateView, lane G 2026-11-13): the state at \a path on the basis these arguments build (as
    //! \c PrepareRestart builds it), then the saved D, C, eps, f -- no Hamiltonian, no SCF.  REFUSED as \c Calculation::JudgeView
    //! refuses, never thrown.  Reports one level in, under the scope its caller opened (a \c report::Section "view").
    static Outcome<std::unique_ptr<StateView>, RestartRefusal>
        ViewState(const std::string& path, int Z, int charge = 0, const RunOptions& opts = {}, const AcceleratorOptions& acc = {});

    bool Converge(const SCFParams& params = {});
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

    size_t           IterationCount() const;
    bool             IsConverged()    const;
    int              GetZ()           const {return itsZ;}
    const Structure& GetStructure()   const {return *itsStructure;}

private:
    typedef std::pair<double, const sf_t*> occ_t;
    struct BuildOnly {};
    //! Everything but the SCF: the members and the electron configuration (the graph is \c BuildGraph's).
    AtomCalculation(BuildOnly, int Z, int charge, const RunOptions& opts, const AcceleratorOptions& acc, Observer onIteration);
    //! The accelerator, the basis (first time), the Hamiltonian and the SCF iterator -- inside the caller's open report run.
    void BuildGraph();
    std::string AcceleratorTag() const;   //!< the accelerator type the caller asked for, as the scf header names it
    bool itsPrepared = false;             //!< \c Prepare built the graph: the next Converge reuses it instead of rebuilding
    void RebuildSampling();
    //! \name THE STATE HOOKS (\c Calculation's template method, V-FACADE-STATE)
    //!@{
    void                    RecordState(StateRecord&, double kT) const override;   //!< the base's keys + the radial basis, model, mesh
    std::string             StateKind() const override {return "atom";}
    const Structure&        StateStructure() const override {return *itsStructure;}
    StateIdentity           MakeStateIdentity() const override;
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

    int                                  itsZ;
    int                                  itsNe;          //!< electron count (Z - charge)
    RunOptions                      itsOpts;
    AcceleratorOptions                   itsAcc;
    std::shared_ptr<const Structure>     itsStructure;   //!< the single atom (Hamiltonian shares it)
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
    std::unique_ptr<qchem::ChargeDensity::rDM_CD> itsCarried;   //!< the converged D the pilot carries into its re-convergence, consumed by \c BuildGraph
    void ResolveGrids();
    //! The electron configuration the SCF occupies: \c opts.configuration STATED (checked, refused when it cannot run) or derived; reported
    //! under policy AtomConfiguration when derived.
    ElectronConfiguration* MakeConfiguration() const;
    void PilotGrids(const SCFParams& params);
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
