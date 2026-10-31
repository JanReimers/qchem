// File: Calculation/Calculation.C
//
// qchem::Calculation -- the production front door for a molecular SCF.  It owns the whole object
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
export module qchem.Calculation;

import qchem.Structure;            // Structure, Molecule, Atom
import qchem.ScalarFunction;       // ::ScalarFunction<double>
import qchem.BasisSet;             // BasisSet::Real_BS
import qchem.Hamiltonian.Factory;  // Hamiltonian::Model, IsDFT, the unified resolver
export import qchem.MoleculeSetup; // CalcOptions + the resolved MoleculeSetup a run is built from (step 7)
export import qchem.PotentialEnergySurface;   // the face MoleculeSession implements (step 8)
export import qchem.SCFFailure;               // SCFFailure: the failure half of the Outcome door (shared with the solid facade)
import qchem.Mesh;                 // qcMesh::MeshParams (the DFT integration grid)
import qchem.ElectronConfiguration;// ElectronConfiguration
import qchem.SCFIterator;          // SCFIterator, SCFParams, IterationRow, EnergyBreakdown
import qchem.Symmetry.Irrep;       // Irrep
import qchem.ChargeDensity;        // rDM_CD
import qchem.ChargeDensity.Seed;   // SeedStrategy
export import qchem.Response;      // the linear response (StaticPolarizability's Outcome and its failure)

export namespace qchem
{

class Calculation
{
public:
    using sf_t       = ScalarFunction<double>;
    using orbitals_t = qchem::Orbitals::Orbitals;
    using Observer   = qchem::SCFIterator::SCFIterator::Observer; //!< void(const IterationRow&): the console's own cells

    //! \brief PROOF that the SCF converged (2026-10-31, the molecular half of the Outcome door -- \c SolidCalculation::Converged's
    //! twin; see its WHY).  \c Converge() and \c Result() hand one out only for a converged run, so a caller that takes its
    //! answers from here cannot read a failed run's number.  The facade's own \c Energy() / \c Density() ... still exist beside it
    //! (the pybind bridge and the GUI read them): moving them onto the proof is the next step, not done here.
    //! \warning NON-OWNING, and invalidated by the next \c Converge() -- take it fresh from each attempt.
    class Converged
    {
    public:
        double                 Energy()         const;   //!< total energy (hartree)
        qchem::EnergyBreakdown EnergyTerms()    const;
        size_t                 IterationCount() const;
        const sf_t&            Density()        const;   //!< rho(r) of the converged state
    private:
        friend class Calculation;
        explicit Converged(const Calculation* c) : itsCalc(c) {}
        const Calculation* itsCalc;
    };

    //! \brief Build the whole graph over a RESOLVED setup (doc/CalcNetworkReview.md §7, step 7) -- the geometry point at the
    //! setup's reference structure -- and converge ONCE on \a scf.  The setup is \c MoleculeSetup::Resolve(structure, options,
    //! accelerator options), which deep-copies the structure (a periodic UnitCell stays periodic), so the caller's may be freed
    //! at once.  The only door (user ruling 2026-10-30: clean replace, no (structure, options) constructors).
    //! \a onIteration, when given, observes THIS constructor's SCF from its first iteration (a GUI streams the trace live);
    //! it stays attached for later \c Converge calls.  A constructor argument, not a setup field: an observer is not a
    //! resolved decision.
    explicit Calculation(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf = {}, Observer onIteration = nullptr);
    ~Calculation();

    Calculation(const Calculation&)            = delete;  // owns raw resources; non-copyable
    Calculation& operator=(const Calculation&) = delete;

    //! Re-run the SCF (e.g. with tighter tolerances) and hand back the ANSWERS or the REASON (it was a bool until 2026-10-31).
    //! Invalidates references previously handed out by Density()/HOMO()/Orbital()/Orbitals() -- the wave function is rebuilt --
    //! and any \c Converged from an earlier attempt.  \c [[nodiscard]] via \c Outcome.
    Outcome<Converged, SCFFailure> Converge(const SCFParams& params = {});
    //! The outcome of the attempt the CONSTRUCTOR made (it converges on construction), so a caller that never re-converges
    //! still faces it.
    Outcome<Converged, SCFFailure> Result() const;

    //! Live per-iteration telemetry.  Set BEFORE a Converge() call to observe it (the ctor's initial
    //! convergence runs before any observer can be attached).
    void OnIteration(Observer obs);

    double                 Energy()      const;   //!< total energy E (hartree)
    qchem::EnergyBreakdown EnergyTerms() const;   //!< the per-term breakdown

    const sf_t&      Density()           const;   //!< rho(r): sample directly, also has Gradient()
    const sf_t&      HOMO()              const;   //!< highest occupied MO
    const sf_t&      Orbital(size_t i)   const;   //!< i-th occupied MO (0 = lowest energy)
    size_t           NumOccupied()       const;
    const orbitals_t* Orbitals(const Irrep&) const; //!< the orbital set for one point-group irrep

    //! One molecular orbital of the converged state: what an orbital panel lists and a field view samples.
    struct Level
    {
        double      e;            //!< eigen-energy (hartree)
        double      occupation;   //!< electrons in it (0 for a virtual)
        Spin        spin;         //!< Up / Down for a polarized run, None for an unpolarized one
        std::string label;        //!< the orbital's QN label (irrep under symmetry)
        const sf_t* phi;          //!< the orbital as a field.  NEVER null; BORROWED from this Calculation
    };
    //! \brief EVERY MO, virtuals included, ascending in energy (both spin channels interleaved for a polarized run).
    //! Cached by each Converge; the vector and every \c Level::phi are valid until the NEXT Converge.
    const std::vector<Level>& GetLevels() const;
    //! \brief m(r) = ρ↑−ρ↓.  An unpolarized run's is the ZERO field (polarized is the general case, unpolarized its
    //! special case), so a caller needs no null branch.  Cached by each Converge; valid until the next.
    const sf_t& GetSpinDensity() const;

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

private:
    //! The finite-difference response ORACLE (src/Response/tests, ruling D6) drives the SAME Reference/frame/
    //! probe/solver as StaticPolarizability with its own kernel, so it needs this run's Hamiltonian, wave
    //! function and basis -- which the facade deliberately does not expose.
    friend class ::ResponseFacadeTests;
    //! \name THE GEOMETRY POINT a MoleculeSession steps to (step 8)
    //!@{
    friend class MoleculeSession;
    //! What crosses from one point to the next: every block's physical D = Σ f c cᵀ, spin irreps outer, basis blocks inner --
    //! the SAME basis recipe at moved atoms has the same blocks, so block b of one point is block b of the next.
    struct Carry { std::vector<hmat_t<double>> D; };
    //! The point at \a at (\c setup->reference's atoms, moved), seeded from \a carry when given (a WorkingDensity), else from the
    //! setup's seed strategy.  The public setup constructor is this at the reference with no carry.
    Calculation(std::shared_ptr<const MoleculeSetup> setup, const Structure& at, std::unique_ptr<const Carry> carry,
                const SCFParams& scf, Observer onIteration);
    Carry MakeCarry() const;
    std::unique_ptr<const Carry> itsCarry;   //!< consumed by the first Converge, which builds the basis it is re-homed onto
    //!@}
    //! The one place a Converged/SCFFailure is minted, so \c Converge() and \c Result() cannot drift.
    Outcome<Converged, SCFFailure> Outcome_() const;
    typedef std::pair<double, const sf_t*> occ_t;   //!< (eigen-energy, orbital) for HOMO/Orbital(i)
    void RebuildSampling();   //!< after a Converge: own a fresh density + sort the occupied MOs

    //! THE RUN SETUP (step 7): the options and every decision resolved from them once.  Everything below is this GEOMETRY POINT
    //! (its structure, EC, basis, Hamiltonian) or its SCF state.
    std::shared_ptr<const MoleculeSetup> itsSetup;
    std::shared_ptr<const Structure>     itsStructure;  //!< the point's structure (valence-charged under a PP); the Hamiltonian shares it
    ElectronConfiguration*               itsEC    = nullptr;  //!< owned
    BasisSet::Real_BS*                   itsBasis = nullptr;  //!< owned
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
    Observer                             itsObserver;          //!< optional live-progress sink
};

//! \brief THE MOLECULAR POTENTIAL ENERGY SURFACE (doc/CalcNetworkReview.md §3.1, §7; step 8): one fixed \c MoleculeSetup,
//! sampled at positions the caller chooses.  Each \c ConvergeAt builds a new point -- the reference's atoms moved, its basis
//! made there (the same named basis) -- seeded from the previous point's density carried onto it, block for block.
//!
//! WHAT A POINT REFUSES (a \c PESFailure, never a crash): positions that change the abelian point group the SALC blocks are
//! built on, when the setup imposes one (\c CalcOptions::symmetry; ruling Q4); forces (step 9).  The molecular basis is
//! neither trimmed nor given an ortho tolerance by this facade, so it has no conditioning window to leave (Q3 is solid-only today).
class MoleculeSession : public virtual PotentialEnergySurface
{
public:
    //! \a onIteration observes every point's SCF (each \c ConvergeAt), from its first iteration.
    MoleculeSession(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf, Calculation::Observer onIteration = nullptr);
    ~MoleculeSession();
    PESGeometry Reference() const override;
    //! THROWS when \a cell is given (a molecule has none) or when \a R has not one position per atom.
    Outcome<PESPoint, PESFailure> ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R,
                                             const PESRequest& want = {}) override;
    //! The current point (the last \c ConvergeAt's); null before the first.  Replaced by the next \c ConvergeAt.
    const Calculation* Current() const {return itsPoint.get();}
private:
    std::shared_ptr<const MoleculeSetup> itsSetup;
    SCFParams                            itsSCF;
    Calculation::Observer                itsObserver;
    std::unique_ptr<Calculation>         itsPoint;
};

} // namespace qchem
