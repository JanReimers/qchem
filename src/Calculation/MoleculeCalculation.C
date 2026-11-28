// File: Calculation/MoleculeCalculation.C (was Calculation.C: renamed 2026-11-02, user ruling -- qchem::Calculation is now the
// abstract base of the three facades)
//
// qchem::MoleculeCalculation -- the production front door for a molecular SCF.  It owns the whole object
// graph (structure copy + EC + basis + Hamiltonian + accelerator + iterator), runs the canonical
// assemble->converge recipe, and exposes only high-level questions: Energy(), Density(), HOMO(),
// Orbital(i).  Density and every MO ARE ScalarFunction<double>s, so a caller samples rho(r) or any
// orbital through one interface (the standout-win that the API-ergonomics review called out).
//
// This is the single tested recipe that previously lived only in the QchemTester test scaffold and
// the pybind bridge's `struct Calc`.  Concrete <double> (molecules) for now; a tCalculation<T>
// extraction is the clean future move when a plane-wave/crystal front door is needed.
module;
#include <memory>
#include <optional>
#include <vector>
#include <string>
#include <utility>
#include "forward.H"   // ResponseFacadeTests -- the FD-kernel response oracle's friend (LinearResponsePlan D6)
export module qchem.MoleculeCalculation;

import qchem.Structure;            // Structure, Molecule, Atom
import qchem.ScalarFunction;       // ::ScalarFunction<double>
import qchem.BasisSet;             // BasisSet::Real_BS
export import qchem.MoleculeSetup; // RunOptions + the resolved MoleculeSetup a run is built from (step 7)
export import qchem.PotentialEnergySurface;   // the face MoleculeSession implements (step 8)
export import qchem.SCFFailure;               // SCFFailure: the failure half of the Outcome door (shared with the solid facade)
export import qchem.Calculation;              // the abstract base of the three facades
export import qchem.StateView;                // ViewState's answer: a saved state, drawn
import qchem.Mesh;                 // qcMesh::MeshParams (the DFT integration grid)
import qchem.ElectronConfiguration;// ElectronConfiguration
import qchem.SCFIterator;          // SCFIterator, SCFParams, IterationRow, EnergyBreakdown
import qchem.Symmetry.Irrep;       // Irrep
import qchem.ChargeDensity;        // rDM_CD
import qchem.ChargeDensity.Seed;   // SeedStrategy
export import qchem.Response.Reference;   // ResponseFailure: StaticPolarizability's failure
export import qchem.LASolver.Krylov;       // KrylovParams (StaticPolarizability)
export import qchem.GridAudit;             // GridAuditResult (MeshAudit)
import qchem.Hamiltonian.Factory;    // CoulombFitSlot (the run's one Coulomb fit, S2b)

export namespace qchem
{

class MoleculeCalculation : public virtual tCalculation<double>
{
public:
    using sf_t       = ScalarFunction<double>;
    using orbitals_t = qchem::Orbitals::Orbitals;
    using Observer   = qchem::SCFIterator::SCFIterator::Observer; //!< void(const IterationRow&): the console's own cells

    //! The proof type is the base's (\c tCalculation<double>::Converged, V-FACADE-PARITY P2): one for every kind.  The facade's own
    //! \c Energy() / \c Density() ... still exist beside it (the pybind bridge and the GUI read them).

    //! NO CONVERGING CONSTRUCTOR (V-FACADE-PARITY R8, user 2026-11-17: "60-90% of the run time should not be in a constructor"): the
    //! one path is Setup (\c MoleculeSetup::Resolve / Admit) -> \c Prepare (the build, an Outcome) -> \c Converge (the Go button).
    ~MoleculeCalculation();

    MoleculeCalculation(const MoleculeCalculation&)            = delete;  // owns raw resources; non-copyable
    MoleculeCalculation& operator=(const MoleculeCalculation&) = delete;

    //! \name BUILT, NOT CONVERGED (2026-11-02: the twin of \c SolidCalculation::Prepare -- the GUI's Build -> Run trigger)
    //! Builds the whole graph over \a setup -- the basis, the Hamiltonian with every mesh it quadratures on, the electron
    //! configuration, the accelerator and the SCF iterator (whose wave-function build files basis.perIrrep / basis.redundant) --
    //! INTO THE CALLER'S report run (a facade never opens one -- user ruling D, 2026-11-02: the GUI's Build opens
    //! "<structure ID> build" itself; with no run open the report writes are no-ops), and stops before iteration 1.  So a client
    //! reads the whole setup BEFORE any SCF.  Run it with \c Converge(params): the FIRST Converge reuses this graph (so the
    //! caller's run around it holds only the SCF), a later one rebuilds as always.  \a onIteration
    //! observes from that first iteration.  Until it runs, \c Result() says "did not converge in 0 iterations" and the energy /
    //! density / orbital accessors have nothing to read.
    //!@{
    //! REFUSED (a value, R9) when a policy refuses while the graph is built (a 4-centre tensor that cannot fit in memory, a PP the run
    //! cannot hold); its policy row is already filed.
    static Outcome<std::unique_ptr<MoleculeCalculation>, BuildRefusal> Prepare(std::shared_ptr<const MoleculeSetup> setup, Observer onIteration = nullptr);
    //!@}

    //! \name RESTART FROM A SAVED STATE (V-FACADE-STATE; \c SaveState / \c JudgeState are the base's)
    //! The molecular twins of \c SolidCalculation::Restart / \c PrepareRestart.  The state's fingerprint is judged against THIS
    //! setup's (\c Calculation::JudgeRestart): REFUSED on another kind, another basis, other electrons or blocks (a value, never a
    //! throw: the caller decides on a fresh run); an EXACT RESUME when nothing differs -- the saved D drives the first Fock
    //! (\c SeedRole::ExactResume) and a Ladder starts on its GDM rung (kT = 0), so the run starts AT the saved answer -- else a WARM
    //! START (\c SeedRole::WorkingDensity), every difference named.  At MOVED atoms the density is carried as AO coefficients
    //! (\c CarryToMovedAtoms, the measured best first iterate).  A setup whose XC mesh is neither stated nor sized ADOPTS the saved state's mesh when the
    //! state was made for the same accuracy eps (its pilot sized it for that eps): no second pilot, and the resume stays exact.
    //!@{
    static Outcome<std::unique_ptr<MoleculeCalculation>, RestartRefusal>
        Restart(const std::string& path, std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf = {}, Observer onIteration = nullptr);
    //! \c Restart up to its SCF (the graph built over the restored density, into the caller's report run); run it with \c Converge.
    static Outcome<std::unique_ptr<MoleculeCalculation>, RestartRefusal>
        PrepareRestart(const std::string& path, std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf = {},
                       Observer onIteration = nullptr);
    //! \c PrepareRestart's judgement, or a \c MoleculeSession point's carry from the previous point (path "(previous geometry
    //! point)", warm, its displacement and the electrons the carried D holds there); empty for a fresh run.
    const std::optional<RestartSummary>& RestartedFrom() const override;
    //!@}

    //! \brief A SAVED STATE, DRAWN (\c StateView, lane G 2026-11-13): the state at \a path on THIS setup's basis -- the basis built
    //! exactly as \c PrepareRestart builds it, then the saved D, C, eps, f; NO Hamiltonian, mesh or fit, no SCF.  REFUSED as
    //! \c Calculation::JudgeView refuses (another kind, a refuse-class difference, ANY moved atom), never thrown.  Writes its
    //! report one level in, under the scope its caller opened (a \c report::Section "view").
    static Outcome<std::unique_ptr<StateView>, RestartRefusal> ViewState(const std::string& path, std::shared_ptr<const MoleculeSetup> setup);

    //! \c Converge(edit) is the base's template method (V-FACADE-PARITY P3): the setup's \c opts.scf, edited by a one-field override
    //! (R10).  It invalidates references previously handed out by Density()/HOMO()/Orbital()/Orbitals() -- the wave function is
    //! rebuilt -- and any \c Converged from an earlier attempt.
    //! \c Result() (the base's): the outcome of the last attempt ("did not converge in 0 iterations" before the first Converge).

    //! Live per-iteration telemetry, REPLACING the one handed to \c Prepare.  Set BEFORE a Converge() call to observe it.
    void OnIteration(Observer obs);

    double                 Energy()      const;   //!< total energy E (hartree)
    qchem::EnergyBreakdown EnergyTerms() const;   //!< the per-term breakdown

    const sf_t&      Density()           const;   //!< rho(r): sample directly, also has Gradient()
    const sf_t&      HOMO()              const;   //!< highest occupied MO
    const sf_t&      Orbital(size_t i)   const;   //!< i-th occupied MO (0 = lowest energy)
    size_t           NumOccupied()       const;
    const orbitals_t* Orbitals(const Irrep&) const; //!< the orbital set for one point-group irrep

    //! One molecular orbital of the converged state: the structure-neutral \c Calculation::Level (e, occupation, spin, k = Γ, label)
    //! PLUS the orbital as a field, which a molecule's field view samples.  Deliberately the same name: a molecule's level IS a level.
    struct Level : Calculation::Level
    {
        const sf_t* phi = nullptr;   //!< the orbital as a field.  NEVER null; BORROWED from this MoleculeCalculation
    };
    //! \brief EVERY MO, virtuals included, in \c MakeLevels' order (ascending in energy, both spin channels interleaved for a polarized
    //! run), so index i here is index i there.  Cached by each Converge; the vector and every \c Level::phi are valid until the NEXT Converge.
    const std::vector<Level>& GetLevels() const;
    //! \brief m(r) = ρ↑−ρ↓.  An unpolarized run's is the ZERO field (polarized is the general case, unpolarized its
    //! special case), so a caller needs no null branch.  Cached by each Converge; valid until the next.
    const sf_t& GetSpinDensity() const;

    //! \name THE CURRENT ITERATE, for a client sampling a run IN PROGRESS (the GUI's "grab 3D plot data", OpenWork §2)
    //! rho(r) / m(r) of the wave function's LATEST diagonalisation (this iteration's output density, before any mixing), BUILT
    //! fresh on each call (an owned copy: one density build, so a client asks only when a person wants a picture).  Call it
    //! from the observer (the compute thread, between iterations: the iterator has updated the wave function before it
    //! reports the row) or between Converge calls -- never from another thread while an SCF runs.  THROWS before the first
    //! iteration (a \c Prepare'd run that has not iterated has no orbitals yet).  Public because the binding samples it from
    //! its row callback, where the cached \c Density() (the CONVERGED state, set after Converge returns) does not exist yet.
    //!@{
    std::unique_ptr<sf_t> MakeIterateDensity()     const;
    std::unique_ptr<sf_t> MakeIterateSpinDensity() const;   //!< the zero field for an unpolarized run
    //!@}

    size_t           IterationCount() const;
    bool             IsConverged()    const;
    //! \brief The static dipole POLARIZABILITY tensor \f$\alpha_{ij}=-\partial^2E/\partial F_i\partial F_j\f$
    //! (bohr³) of the converged state, by a SELF-CONSISTENT linear response on the Hamiltonian's analytic kernel
    //! -- CPHF for HF (doc/LinearResponsePlan.md stage R1).  The dipole matrices are numerical, on
    //! \a dipoleMesh (ruling Q5).  MEASURED on H2O/dzvp: the numerical dipole oscillates about the analytic
    //! answer at ~1e-7 relative whatever the mesh (a Becke-quadrature floor; IntegrationTests/M_Response.C).
    //! FAILS on an ungapped reference (E1) or a Krylov solve that does not reach \a krylov.tol.  THROWS for a
    //! model whose Hamiltonian cannot be linearised yet (LDA until R2 gives its fitted terms the face).
    //! Reports the per-channel χ0 / χ and residuals on stdout, at its own activity.
    Outcome<rmat_t,Response::ResponseFailure> StaticPolarizability(
        const qcMesh::MeshParams& dipoleMesh = {.radial=qcMesh::RadialKind::MHL, .nRadial=80, .mhl_m=3, .mhl_alpha=2.0,
                                                .angular=qcMesh::AngularKind::Lebedev, .angularDegree=35, .beckeOrder=3},
        const KrylovParams& krylov = {.tol=1e-10}) const;
    const Structure& GetStructure()   const {return *itsStructure;}
    //! \brief THE ENERGY-WEIGHTED DENSITY W = Σ f ε c cᵀ of the converged point, per basis block, in the carry's order (spin irreps outer,
    //! basis blocks inner), so Σ_blocks Tr(W S) = Σ_levels f ε.  PUBLIC, and why (CLAUDE.md: say why at the declaration): its client is the
    //! analytic force's orthonormality (Pulay) term +Tr(W ∂S/∂R), which the force assembly needs above the wave function
    //! (doc/GeometryOptimizationPlan.md §13: W is the facade's, WaveFunction gains nothing); a cached Get would hold matrices no other
    //! client wants.  THROWS before a Converge (no wave function).
    std::vector<hmat_t<double>> MakeEnergyWeightedDensity() const;
    //! \brief THE XC-MESH AUDIT (doc/GeometryOptimizationPlan.md §14.3; pin 31): the run against a Hamiltonian rebuilt on the SAME basis
    //! with every mesh-borne term on \a reference (the XC mesh, and a δ-ρ Hartree rung that rides it), at the converged D, term by term.
    //! A molecule's basis does not depend on its mesh, so nothing is restored: the one density serves both sides.  Terms the mesh does not
    //! touch (kinetic, Ven, the GTH terms, Dunlap J, ion-ion) read 0.  The run's own Hamiltonian is re-refreshed for that D.  Writes the
    //! console table and its fields AT THE REPORT CURSOR: the caller opens the scope (the report nesting rule); the reference
    //! Hamiltonian's own mesh reports go to a nested run, so this run's `hamiltonian` keeps describing the meshes it USED.  Judged on
    //! \c Consumes(): |ΔE| against eps x N_atoms, or max|ΔH_ij| against eps.  THROWS before a Converge.
    GridAuditResult MeshAudit(const qcMesh::MeshParams& reference) const;
    //! \brief THE FOLD AUDIT (the per-orbit consumer of the molecular mesh fold, OpenWork §3): the run against a Hamiltonian rebuilt on the
    //! SAME mesh with the OPPOSITE fold choice (folded <-> unfolded), at the converged D, term by term, as \c MeshAudit.  The fold is exact,
    //! not an approximation, so every term reads roundoff: the gate that the per-orbit route computes the same Hamiltonian cheaper.
    //! THROWS (invalid_argument) on a free run (its density is not invariant: the folded mesh is a different answer) or a run with no
    //! point group on its mesh; THROWS before a Converge.
    GridAuditResult FoldAudit() const;
    //! \brief One Hartree route against another at ONE density (\c CoulombAudit).  \c ΔE_H = E_H(run) - E_H(reference); \c maxAbsΔJ =
    //! max over every block, spin and element of |F(run) - F(reference)|, which is |J(run) - J(reference)| because nothing else differs.
    struct CoulombAuditResult
    {
        std::string reference;           //!< the reference route, spelled as the deck states it
        double      EH_run=0.0, EH_ref=0.0, ΔE_H=0.0, maxAbsΔJ=0.0;
        double      wallSeconds=0.0;     //!< the reference's build and one Fock, at this D
    };
    //! \brief THE COULOMB AUDIT (lane S 2026-11-25, the Q2 fixed-D column; controller ruling the same day): the run's Hartree against
    //! \a reference, at the converged D, no SCF.  The reference Hamiltonian is built on the run's own mesh and mesh symmetry with only the
    //! Coulomb route changed, so every other term is the same term and the Fock difference is the J difference; the audit CHECKS that
    //! (every non-Hartree energy equal to 1e-10, else logic_error).  \a reference = \c HartreeRoute::Exact is the accuracy oracle (the
    //! four-centre J).  An \c Auxiliary reference must state its \c auxBasis (no policy runs here).  Writes its fields AT THE REPORT CURSOR
    //! (the caller opens the scope); the reference's own reports go to a nested run.  THROWS before a Converge, and (invalid_argument) on a
    //! route \c AdmitCoulomb refuses.
    CoulombAuditResult CoulombAudit(const qchem::Hamiltonian::CoulombRoute& reference) const;
    //! \brief WHAT THIS RUN'S CLIENT CONSUMES (pin 31; the solid's rule): a run built on its own answers an \c Energy, judged on |ΔE| against
    //! eps x N_atoms; a point a \c MoleculeSession builds (a surface, forces, an optimization) answers \c MatrixElements, max|ΔH_ij| against eps.
    ClientConsumes Consumes() const override;
    //! \brief THE MESH PILOT's cap on climb + re-convergence rounds (the solid's \c kPilotRounds): round 0 climbs at the first converged D,
    //! each later round is the CONFIRMING audit (one raise at most) at the re-converged D.
    static constexpr int kPilotRounds=4;
    //! \brief What the mesh pilot COST, reported, never hidden (doc/GeometryOptimizationPlan.md §14.2 step 5): its rounds, the reference
    //! Hamiltonians it built, and the SCF iterations of its re-convergences (on top of the run's own).  settled: the climb met the target
    //! on both axes (not the round cap, not a ladder end, not a stated mesh over it).  statedOver: a STATED mesh missed the target --
    //! honoured, reported (a `warn` policy row), and not settled (pin 31 rule 3).  All zero on a run with no mesh (HF) or a resolved one.
    struct MeshPilotCost { int rounds=0; int referencesBuilt=0; size_t reconvergeIterations=0; bool settled=false; bool statedOver=false; };
    MeshPilotCost PilotCost() const {return itsPilotCost;}
    //! \brief The setup with THIS run's XC mesh FROZEN into it: the rung the pilot settled on, marked resolved -- what a \c MoleculeSession
    //! steps the rest of its points on, so every point integrates on the same mesh (a re-sized mesh per point would make the surface
    //! discontinuous; §2 rule 4).  A stated mesh stays stated.  The setup itself when nothing was piloted.
    std::shared_ptr<const MoleculeSetup> SizedSetup() const {return itsSetup;}

private:
    //! The finite-difference response ORACLE (src/Response/tests, ruling D6) drives the SAME Reference/frame/
    //! probe/solver as StaticPolarizability with its own kernel, so it needs this run's Hamiltonian, wave
    //! function and basis -- which the facade deliberately does not expose.
    friend class ::ResponseFacadeTests;
    //! The per-term finite difference at FIXED D (doc/GeometryOptimizationPlan.md §6): a term's energy at moved atoms with the carried
    //! density frozen, which needs this point's Hamiltonian and the carried density on its basis -- both private.
    friend class ::ForceFacadeTests;
    friend class ::GeometryCacheTests;   // S2b's gate reads a point's Coulomb fit and basis blocks
    //! \name THE GEOMETRY POINT a MoleculeSession steps to (step 8)
    //!@{
    friend class MoleculeSession;
    //! What crosses from one point to the next: every block's physical D = Σ f c cᵀ, spin irreps outer, basis blocks inner --
    //! the SAME basis recipe at moved atoms has the same blocks, so block b of one point is block b of the next.
    //! A SESSION STEP fills \c moved instead (\c D empty): the previous point's state, carried onto this point's basis by \c BuildGraph
    //! once that basis exists, through the SAME function a file restart at moved atoms uses (\c CarryToMovedAtoms), which files the
    //! step in \c RestartedFrom.  \c D is the plain AO carry, for the same geometry (the pilot's re-convergence) and for the fixed-D
    //! probes of a finite difference, which DEFINE the frozen D as an AO matrix.
    struct Carry
    {
        std::vector<hmat_t<double>> D;
        std::optional<SavedState>   moved;              //!< a session step: D and S of the previous point (\c MakeMovedCarry)
        double                      displacement = 0.0; //!< the step's largest atom move (bohr), for \c RestartedFrom
    };
    //! The point \a point (\c MoleculeSetup::MakePoint: the setup's reference atoms, moved, with the mesh maker over them), seeded from
    //! \a carry when given (a WorkingDensity), else from the setup's seed strategy.  \c Prepare is this at the reference with no carry; a
    //! \c MoleculeSession point, this at its atoms.
    struct BuildOnly {};
    //! Everything but the SCF: the members and the electron configuration (the graph is \c BuildGraph's).
    MoleculeCalculation(BuildOnly, std::shared_ptr<const MoleculeSetup> setup, MoleculeSetup::Point point, std::unique_ptr<const Carry> carry,
                Observer onIteration, bool geometryPoint=false);
    //! The accelerator, the basis (first time), the Hamiltonian and the SCF iterator -- inside the caller's open report run.
    void BuildGraph();
    bool itsPrepared = false;   //!< \c Prepare built the graph: the next Converge reuses it instead of rebuilding
    //! The point's Hamiltonian on \a mesh (the setup's model or pseudopotential, this basis): the graph's, and every audit reference's.
    qchem::Hamiltonian::rHamiltonian* MakeHamiltonianOn(const qcMesh::MeshParams& mesh) const;
    //! The same with the mesh source \a symmetry in place of the point's (the fold audit's reference).
    qchem::Hamiltonian::rHamiltonian* MakeHamiltonianOn(const qcMesh::MeshParams& mesh, const qcMesh::SymmetrySource& symmetry) const;
    //! The same with the Hartree route \a coulomb and its fit slot \a fit in place of the setup's (the Coulomb-fit check's δ reference), and
    //! the mesh source \a symmetry when given (null: the point's).
    qchem::Hamiltonian::rHamiltonian* MakeHamiltonianOn(const qcMesh::MeshParams& mesh, const qchem::Hamiltonian::CoulombRoute& coulomb,
                                                        const qchem::Hamiltonian::CoulombFitSlot<double>& fit,
                                                        const qcMesh::SymmetrySource* symmetry=nullptr) const;
    //! The audit body: the run against the Hamiltonian on \a reference built under \a symmetry, named \a referenceName, under \a title.
    GridAuditResult AuditAgainst(const qcMesh::MeshParams& reference, const qcMesh::SymmetrySource& symmetry, const std::string& referenceName,
                                 const char* title) const;
    //! THE COULOMB-FIT CHECK (policy CoulombFitResidual): the fit's J error, J_δ − J_fit, at \a d -- the seed's first density before the
    //! SCF (\a when "seed"), the converged one after (\a when "converged").  WARNs over the accuracy target; skipped, saying so, when the
    //! δ reference would not fit in memory.
    void CheckCoulombFit(const qchem::ChargeDensity::rDM_CD* d, const std::string& when);
    std::string FitCaveat(double target) const;   //!< the mesh pilot's "meets" line on a Gaussian-fit Hartree: the fit's error is not the mesh's
    //! Does this run integrate on a mesh at all (a DFT model or a pseudopotential)?  HF has none: nothing to audit or climb.
    bool HasXCMesh() const;
    //! THE MESH PILOT (§14.2): climb at the converged D, rebuild + re-converge once per raise, freeze into \c itsSetup.  From \c Converge.
    void PilotMesh(const SCFParams& params);
    //! Point the new iterator's observer at this facade (the last row, and the caller's sink).
    void AttachObserver();
    bool          itsGeometryPoint = false;   //!< built by a MoleculeSession: answers MatrixElements (\c Consumes)
    bool          itsPiloted       = false;   //!< the pilot ran (once per calculation: a later Converge keeps its mesh)
    bool          itsFitConfirmed  = false;   //!< the Coulomb-fit check ran on a converged density (once per calculation)
    double        itsFitResidual   = -1.0;    //!< the latest Coulomb-fit check's J error (Ha); < 0 = not measured.  The mesh pilot's line cites it
    MeshPilotCost itsPilotCost;
    size_t itsIterations=0;   //!< the last Converge's SCF iterations, the pilot's re-convergences included (the iterator counts only its own)
    Carry MakeCarry() const;
    //! The carry a session hands the NEXT point: this point's blocks (\c SnapshotState), carried there by \c CarryOntoThisBasis.
    Carry MakeMovedCarry(double displacement) const;
    //! \a carry's \c moved state taken onto THIS point's basis (\c CarryToMovedAtoms), as an AO carry; files what it did in
    //! \c itsRestartedFrom (a warm start from the previous point, with the electrons it holds here).
    Carry CarryOntoThisBasis(const Carry& carry);
    //! Σ f c cᵀ (or Σ f ε c cᵀ, \a energyWeighted) per block, in the carry's order: one loop for the carry and for W.
    std::vector<hmat_t<double>> BlockDensities(bool energyWeighted) const;
    //! \a carry's blocks as a density on THIS point's basis (throws when the blocks do not match): the carry's WorkingDensity, and the
    //! frozen D of a fixed-D finite difference.
    std::unique_ptr<qchem::ChargeDensity::rDM_CD> MakeCarriedDensity(const Carry& carry) const;
    //! THE POINT'S FORCE (doc/GeometryOptimizationPlan.md §13), −dE/dR in Ha/bohr in the structure's atom order: every term's force at
    //! the converged D plus the Pulay sum Tr(W ∂S/∂R) at its energy-weighted density W.  For the session's forces request.  THROWS (the
    //! Hamiltonian's std::logic_error) when \c MissingGradients is not empty, and before a Converge.
    rvec3vec_t MakeForces() const;
    //! The terms (and the Pulay sum) that cannot give their geometry gradient: the BUILT Hamiltonian's answer, no density needed.
    std::vector<std::string> MissingGradients() const;
    std::unique_ptr<const Carry> itsCarry;  //!< consumed by the first Converge, which builds the basis it is re-homed onto
    //!@}
    //! \name THE STATE HOOKS (\c Calculation's template method, V-FACADE-STATE)
    //!@{
    void                    RecordState(StateRecord&, double kT) const override;   //!< the base's keys + the basis, PP, model, mesh
    std::string             StateKind() const override {return "molecule";}
    const Structure&        StateStructure() const override {return *itsStructure;}
    std::string             StateFunctional() const override;
    std::string             KindTag() const override {return "molecule";}
    double                  StateKT() const override;
    std::vector<StateBlock> StateBlocks(bool withOrbitals) const override;
    StateSummary            MakeStateSummary() const override;
    //!@}
    std::optional<RestartSummary> itsRestartedFrom;   //!< what PrepareRestart (or a session's moved carry) did, empty for a fresh run
    bool   itsCarryExact = false;   //!< \c itsCarry is an EXACT RESUME (\c PrepareRestart): the first Fock is built from it
    bool   itsStartOnGDM = false;   //!< the next accelerator built is a Ladder started on its GDM rung (an exact resume at kT = 0), once
    double itsLastkT     = 0.0;     //!< the last Converge's smearing (the saved state's scf.kT)
    //! \name \c Converge's STEPS (\c tCalculation's hooks)
    //!@{
    const Setup& GetSetup() const override;
    //! One stage: the recipe's report, the graph (\c Prepare's, reused once, else rebuilt), the SCF, the mesh pilot, the result line.
    bool RunStage(const SCFParams& params) override;
    //! A schedule's stage (P4): on the setup's accelerator and recipe (\c AdmitStageOnSetup refuses another), from the density the last
    //! stage reached once any SCF has run -- stage 0 of a \c Prepare'd graph runs on that graph, from the setup's seed.
    void EnterStage(const SCFStage& stage, size_t index) override;
    //!@}
    //! \name THE PROOF'S ANSWERS (\c tCalculation's hooks): the one place a failure is ranked, so \c Converge() and \c Result() cannot drift
    //!@{
    std::optional<SCFFailure>             MakeFailure()      const override;
    EnergyBreakdown                       ProofEnergyTerms() const override;
    double                                ProofCharge()      const override;
    size_t                                ProofIterations()  const override;
    const ChargeDensity::tDM_CD<double>&  ProofDensity()     const override;
    const ScalarFunction<double>&         ProofSpinDensity() const override;
    //!@}
    typedef std::pair<double, const sf_t*> occ_t;   //!< (eigen-energy, orbital) for HOMO/Orbital(i)
    void RebuildSampling();   //!< after a Converge: own a fresh density + sort the occupied MOs

    //! THE RUN SETUP (step 7): the options and every decision resolved from them once.  Everything below is this GEOMETRY POINT
    //! (its structure, EC, basis, Hamiltonian) or its SCF state.
    std::shared_ptr<const MoleculeSetup> itsSetup;
    std::shared_ptr<const Structure>     itsStructure;  //!< the point's structure (valence-charged under a PP); the Hamiltonian shares it
    std::shared_ptr<const qcMesh::MeshMaker> itsMeshMaker;  //!< builds meshes over \c itsStructure (StructureMeshMaker of its concrete type)
    std::shared_ptr<const qcMesh::SymmetrySource> itsMeshSymmetry;  //!< what \c itsStructure's meshes respect (the point's, policy MeshSymmetry)
    ElectronConfiguration*               itsEC    = nullptr;  //!< owned
    BasisSet::Real_BS*                   itsBasis = nullptr;  //!< owned
    //! The run's ONE Coulomb (Dunlap) fit on \c itsBasis (doc/IntegralCacheOwnershipPlan.md S2b): the Hartree term of the first Hamiltonian
    //! fills it and every later one (the pilot ladder's rungs, the final build, an audit reference) reads it, so they share its metric
    //! and inverse.  It owns no mesh, and this point's structure and basis never change, so it is a pure function of fixed inputs.
    qchem::Hamiltonian::CoulombFitSlot<double> itsCoulombFit;
    //! A SALC run's UNSYMMETRIZED basis (the adapted one keeps it alive too; null without symmetry): its shells identify a saved
    //! state (RecordState), which a SALC block, a combination of them, cannot report.
    std::shared_ptr<const BasisSet::Real_BS> itsAoBasis;
    //! R2.22: the iterator no longer deletes what it is handed, so this facade owns the Hamiltonian and
    //! the accelerator.  Both are still REBUILT per Converge() -- a molecular Hamiltonian is cheap, so
    //! that stays as it was; assigning here frees the previous pair at exactly the old moment.  Declared
    //! BEFORE itsScf so they outlive the iterator built over them.
    std::unique_ptr<qchem::Hamiltonian::rHamiltonian>       itsHam;
    std::unique_ptr<qchem::SCFAccelerators::SCFAccelerator> itsAccel;
    qchem::SCFIterator::SCFIterator*     itsScf   = nullptr;  //!< owned (the Hamiltonian/accelerator are NOT its)
    std::unique_ptr<qchem::ChargeDensity::rDM_CD> itsDensity;  //!< owned converged rho(r)
    std::unique_ptr<const sf_t>          itsSpinDensity;       //!< owned converged m(r), the zero field when unpolarized
    std::vector<Level>                   itsLevels;            //!< every MO, ascending; phi into the wave function
    std::vector<occ_t>                   itsOccupied;          //!< non-owning, into the wave function
    Observer                             itsObserver;          //!< optional live-progress sink (the last iterate's row is the base's: KeepingLastRow)
};

//! \brief THE MOLECULAR POTENTIAL ENERGY SURFACE (doc/CalcNetworkReview.md §3.1, §7; step 8): one fixed \c MoleculeSetup,
//! sampled at positions the caller chooses.  Each \c ConvergeAt builds a new point -- the reference's atoms moved, its basis
//! made there (the same named basis) -- seeded from the previous point's density carried onto it, block for block.
//!
//! WHAT A POINT REFUSES (a \c PESFailure, never a crash): positions that change the abelian point group the SALC blocks are
//! built on, when the setup imposes one (\c RunOptions::symmetry; ruling Q4); forces while the point's Hamiltonian names a term
//! without its geometry gradient (\c Hamiltonian::MissingGradients: a GGA's XC term, or the δ-mesh Hartree until its operator pullback lands), refused BEFORE any SCF; stress (a molecule has
//! no cell).  A forces request otherwise carries −dE/dR (the terms' forces and the Pulay sum at the converged D and W) at every point.
//! The molecular basis is neither trimmed nor given an ortho tolerance by this facade, so it has no conditioning window to leave (Q3 is
//! solid-only today).
class MoleculeSession : public virtual PotentialEnergySurface
{
public:
    //! \a onIteration observes every point's SCF (each \c ConvergeAt), from its first iteration.
    MoleculeSession(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf, MoleculeCalculation::Observer onIteration = nullptr);
    ~MoleculeSession();
    PESGeometry Reference() const override;
    //! THROWS when \a cell is given (a molecule has none) or when \a R has not one position per atom.
    Outcome<PESPoint, PESFailure> ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R,
                                             const PESRequest& want = {}) override;
    //! The current point (the last \c ConvergeAt's); null before the first.  Replaced by the next \c ConvergeAt.
    const MoleculeCalculation* Current() const {return itsPoint.get();}
private:
    std::shared_ptr<const MoleculeSetup> itsSetup;
    SCFParams                            itsSCF;
    MoleculeCalculation::Observer                itsObserver;
    std::unique_ptr<MoleculeCalculation>         itsPoint;
    size_t                                       itsPointsRun=0;   //!< points this session converged: the timing breakdown's point index
};

} // namespace qchem
