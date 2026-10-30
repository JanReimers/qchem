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
import qchem.Mesh;                 // qcMesh::MeshParams (the DFT integration grid)
import qchem.ElectronConfiguration;// ElectronConfiguration
import qchem.SCFIterator;          // SCFIterator, SCFParams, SCFProgress, EnergyBreakdown
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
    using Observer   = qchem::SCFIterator::SCFIterator::Observer; //!< void(const SCFProgress&)

    //! Build the whole graph and converge.  \a st is DEEP-COPIED via Structure::Clone (which PRESERVES the
    //! concrete geometry -- a periodic UnitCell stays periodic, not sliced to a Molecule), so the caller's
    //! Molecule/Atom/UnitCell is left untouched and may be freed immediately.
    explicit Calculation(const Structure&          st,
                         const CalcOptions&         opts = {},
                         const AcceleratorOptions&  acc  = {});
    //! Build the whole graph and converge ONCE, on \a scf (CleanCode D-CALC-DOUBLE-SCF): a caller with its own recipe used to
    //! construct (an SCF on the defaults) and then Converge(par) -- two full SCFs, and the first one's cap was out of its reach.
    Calculation(const Structure&          st,
                const CalcOptions&         opts,
                const SCFParams&           scf,
                const AcceleratorOptions&  acc  = {});
    //! \brief A run over a RESOLVED setup (doc/CalcNetworkReview.md §7, step 7): the geometry point at the setup's reference
    //! structure, converged once on \a scf.  The constructors above are this with \c MoleculeSetup::Resolve in front.
    explicit Calculation(std::shared_ptr<const MoleculeSetup> setup, const SCFParams& scf = {});
    ~Calculation();

    Calculation(const Calculation&)            = delete;  // owns raw resources; non-copyable
    Calculation& operator=(const Calculation&) = delete;

    //! Re-run the SCF (e.g. with tighter tolerances).  Invalidates references previously handed out
    //! by Density()/HOMO()/Orbital()/Orbitals() -- the wave function is rebuilt.
    bool Converge(const SCFParams& params = {});

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
    std::vector<occ_t>                   itsOccupied;          //!< non-owning, into the wave function
    Observer                             itsObserver;          //!< optional live-progress sink
};

} // namespace qchem
