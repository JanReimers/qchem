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
export module qchem.Calculation;

import qchem.Structure;            // Structure, Molecule, Atom
import qchem.ScalarFunction;       // ::ScalarFunction<double>
import qchem.BasisSet;             // BasisSet::Real_BS
import qchem.Hamiltonian.Factory;  // Hamiltonian::Model, Hamiltonian::Pol, IsDFT, the unified resolver
import qchem.Mesh;                 // qcMesh::MeshParams (the DFT integration grid)
import qchem.ElectronConfiguration;// ElectronConfiguration
import qchem.SCFIterator;          // SCFIterator, SCFParams, SCFProgress, EnergyBreakdown
import qchem.Symmetry.Irrep;       // Irrep
import qchem.ChargeDensity;        // rDM_CD
import qchem.ChargeDensity.Seed;   // SeedStrategy

export namespace qchem
{

using Hamiltonian::Model;          // {E1, HF, DE1, DHF}
using Hamiltonian::Pol;            // {UnPolarized, Polarized}

//! Orbital-integral engine: the in-house MnD recursion (default) or the libcint foreign engine.
enum class Engine  { MnD, LibCint };
//! Angular basis convention: Cartesian Gaussians (default) or real solid-harmonic (spherical).
enum class Angular { Cartesian, Spherical };

//! How to set up the calculation.  Designated-initializer friendly:
//!     Calculation calc(water, {.basis="dzvp"});
//! `model` selects the Hamiltonian; pass 1 wires the non-DFT factory path (HF is the default).
struct CalcOptions
{
    std::string basis = "sto-3g";
    Model       model = Model::HF;   //!< HF (default) | Xalpha | LDA | E1/DE1/DHF (test-only)
    Pol         pol   = Pol::UnPolarized;
    //! Spin multiplicity 2S+1.  0 (default) = minimal spin: closed-shell singlet for even Ne, doublet for
    //! odd -- the historical behaviour.  Set explicitly for an open shell: 3 = triplet, 2 = doublet, ...
    //! The facade converts it to (nUp,nDown) [nUp-nDown = 2S = multiplicity-1, nUp+nDown = Ne] and PROMOTES
    //! the calculation to Pol::Polarized when 2S>0 (unrestricted open shell needs distinct up/down densities).
    //! A multiplicity whose parity disagrees with Ne (e.g. a singlet for odd Ne) is rejected.
    int         multiplicity = 0;
    //! Basis construction variants (threaded into BasisSet::Molecule::Factory).  Defaults reproduce
    //! today's behaviour (in-house MnD, Cartesian).  angular==Spherical + symmetry is rejected until
    //! the Spherical SALC track (doc/SphericalSALCPlan.md) lands -- the SALC builder needs Cartesian PGData.
    Engine      engine  = Engine::MnD;
    Angular     angular = Angular::Cartesian;
    //! DFT-only knobs (ignored when model is HF/1-e/Dirac).  xalpha: the Slater exchange parameter, used
    //! only by model==Xalpha.  mesh: the numerical XC integration grid -- defaults to the proven molecular
    //! values; a designated initializer overrides just the resolution you care about, e.g. {.nRadial=50}.
    double      xalpha = 0.7;
    qcMesh::MeshParams mesh = {.radial  = qcMesh::RadialKind::MHL,   .nRadial   = 30,
                              .mhl_m    = 3,                         .mhl_alpha = 2.0,
                              .angular  = qcMesh::AngularKind::Gauss, .nAngular  = 12,
                              .beckeOrder = 2};
    //! Point-group SALC blocking + per-irrep aufbau.  GUARDED TO THE CARTESIAN PG BASIS: the SALC
    //! builder needs a PolarizedGaussian (PGData) orbital IBS and throws otherwise.  Since the facade
    //! only builds the default Cartesian basis today, this is always the supported path; the guard
    //! future-proofs the day spherical/libcint deliveries are exposed (see doc/SphericalSALCPlan.md).
    bool        symmetry    = false;
    double      symmetryTol = 1e-4;   //!< geometry tolerance for point-group detection
    //! Replace the bare nuclear attraction with the GTH/HGH pseudopotential (valence-only): the atoms keep
    //! only their Zion valence electrons and the Hamiltonian gets V_loc(r) + the KB-separable projectors +
    //! Zion ion-ion, in place of Ven.  SINGLE-SPECIES for now (all atoms the same element; multi-species is
    //! the next increment).  Requires a valence basis in \c basis (an all-electron basis has spurious core
    //! shells under a PP) -- e.g. {.basis="sipp", .pseudopotential=true} for silicon.  Ignores \c model
    //! (the PP front door is always LSDA today).  See doc/MolecularPseudopotentialPlan.md.
    bool        pseudopotential = false;
    //! SCF seed.  Default == auto: DFT seeds from SAD (superposition of atomic densities), HF/1-e take
    //! the core guess.  Set explicitly to override -- e.g. {.seed=SAD} drives HF's matrix-free-seed
    //! bootstrap (the iterator manufactures a D0 via a one-step LDA sibling).
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Default;
};

//! DIIS accelerator knobs.  These were stringly-typed json keys grepped out of tests; here they
//! are documented fields with defaults.  eMax<=0 means "auto" -- the Z-scaled heuristic the
//! molecular recipe has always used.  The long tail still rides through the json factory directly.
struct AcceleratorOptions
{
    int    nProj = 4;
    double eMax  = 0.0;     //!< <=0 => Z*Z*0.1/32 (the proven molecular default)
    double eMin  = 1e-7;
    double svTol = 5e-9;
    //! SCF accelerator: "DIIS" (default) | "GDM" (direct minimisation, robust for hard/open-shell SCF) |
    //! "Ladder" (auto None->DIIS->GDM).  A hard pseudopotential (open-shell atom, diffuse valence) that
    //! limit-cycles under DIIS often converges under GDM or Ladder.
    std::string type = "DIIS";
};

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
    const Structure& GetStructure()   const {return *itsStructure;}

private:
    typedef std::pair<double, const sf_t*> occ_t;   //!< (eigen-energy, orbital) for HOMO/Orbital(i)
    void RebuildSampling();   //!< after a Converge: own a fresh density + sort the occupied MOs

    std::shared_ptr<const Structure>     itsStructure;  //!< facade's own copy (Hamiltonian shares it)
    CalcOptions                          itsOpts;
    AcceleratorOptions                   itsAcc;
    ElectronConfiguration*               itsEC    = nullptr;  //!< owned
    BasisSet::Real_BS*                   itsBasis = nullptr;  //!< owned
    qchem::SCFIterator::SCFIterator*     itsScf   = nullptr;  //!< owned (owns Hamiltonian + accelerator)
    std::unique_ptr<qchem::ChargeDensity::rDM_CD> itsDensity;  //!< owned converged rho(r)
    std::vector<occ_t>                   itsOccupied;          //!< non-owning, into the wave function
    Observer                             itsObserver;          //!< optional live-progress sink
};

} // namespace qchem
