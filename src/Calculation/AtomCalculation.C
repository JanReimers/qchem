// File: Calculation/AtomCalculation.C
//
// qchem::AtomCalculation -- the production front door for a single-atom SCF (the atom/Dirac/PP analogue
// of qchem::Calculation, which is molecule-only by deliberate design).  It owns the whole object graph
// (atom + atomic EC + atomic basis + Hamiltonian + accelerator + iterator), runs the assemble->converge
// recipe, and exposes high-level questions: Energy(), EnergyTerms(), Density(), Orbital(i).
//
// This is the single tested atom recipe that previously lived only in the QchemTester/TestAtom test
// scaffold.  A SIBLING of qchem::Calculation rather than an extension: the molecular facade is being
// consumed by the GUI team and must stay stable, and atoms genuinely differ (atomic exponent-pool bases,
// atomic / Dirac / pseudo ECs, the nAngular=1 atomic mesh, relativistic + PP models the molecule facade
// has no business expressing).  The two share the SCFIterator converge machinery and AcceleratorOptions.
module;
#include <memory>
#include <vector>
#include <string>
#include <utility>
export module qchem.AtomCalculation;

import qchem.Structure;             // Structure, Atom
import qchem.ScalarFunction;        // ScalarFunction<double>
import qchem.BasisSet;              // BasisSet::Real_BS
import qchem.BasisSet.Atom.Factory; // BasisSet::Atom::Type, BasisSetAccuracy (the atomic exponent-pool basis)
import qchem.Hamiltonian.Factory;   // Hamiltonian::Model, Hamiltonian::Pol
import qchem.Mesh;                  // qcMesh::MeshParams
import qchem.ElectronConfiguration; // ElectronConfiguration
import qchem.SCFIterator;           // SCFIterator, SCFParams, SCFProgress, EnergyBreakdown
import qchem.Symmetry.Irrep;        // Irrep
import qchem.Symmetry.Spin;         // Spin (per-spin irrep enumeration)
import qchem.Orbitals;              // Orbital, Orbitals
import qchem.ChargeDensity;         // DM_CD
import qchem.ChargeDensity.Seed;    // SeedStrategy
import qchem.Calculation;           // reuse AcceleratorOptions (read-only; does NOT modify Calculation)

export namespace qchem
{

using Hamiltonian::Model;                 // {E1, HF, DE1, DHF, Xalpha, LDA}
using Hamiltonian::Pol;                   // {UnPolarized, Polarized}
using AtomType = BasisSet::Atom::Type;    // {Slater, Gaussian, BSpline6, BSpliner6, Gaussian_RKB, Slater_RKB}
using BasisSet::Atom::BasisSetAccuracy;   // {N3, N5, Low, Medium, High}

//! How to set up a single-atom calculation.  Designated-initializer friendly:
//!     AtomCalculation calc(18, 0, {.type=AtomType::Slater, .accuracy=BasisSetAccuracy::Medium});
struct AtomCalcOptions
{
    AtomType         type     = AtomType::Slater;            //!< atomic basis family
    BasisSetAccuracy accuracy = BasisSetAccuracy::Medium;    //!< preset exponent pool (used when N<=0)
    //! Explicit exponent pool -- when N>0 these override the accuracy preset (the {type,N,emin,emax} json
    //! the scaffold built by hand for the A_SG/A_SL exponent-sweep tests).
    int    N    = 0;
    double emin = 0.0;
    double emax = 0.0;

    Model  model  = Model::HF;          //!< HF (default) | E1 | DE1/DHF (Dirac) | Xalpha | LDA
    Pol    pol    = Pol::UnPolarized;
    double xalpha = 0.7;                 //!< Slater-Xalpha exchange parameter (model==Xalpha only)

    //! Atomic XC integration grid: the proven atom values (one angular point -- atoms are spherical).
    qcMesh::MeshParams mesh = {.radial  = qcMesh::RadialKind::MHL,   .nRadial   = 50,
                              .mhl_m    = 3,                         .mhl_alpha = 2.0,
                              .angular  = qcMesh::AngularKind::Gauss, .nAngular  = 1,
                              .beckeOrder = 2};
    //! SCF seed.  Default == auto: CoreGuess (atoms never use the molecular SAD seed).
    qchem::ChargeDensity::SeedStrategy seed = qchem::ChargeDensity::SeedStrategy::Default;
};

class AtomCalculation
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
                             const AtomCalcOptions&    opts   = {},
                             const SCFParams&          params = {},
                             const AcceleratorOptions& acc    = {});
    ~AtomCalculation();

    AtomCalculation(const AtomCalculation&)            = delete;
    AtomCalculation& operator=(const AtomCalculation&) = delete;

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
    void RebuildSampling();

    int                                  itsZ;
    int                                  itsNe;          //!< electron count (Z - charge)
    AtomCalcOptions                      itsOpts;
    AcceleratorOptions                   itsAcc;
    std::shared_ptr<const Structure>     itsStructure;   //!< the single atom (Hamiltonian shares it)
    ElectronConfiguration*               itsEC    = nullptr;  //!< owned
    BasisSet::Real_BS*                   itsBasis = nullptr;  //!< owned
    qchem::SCFIterator::SCFIterator*     itsScf   = nullptr;  //!< owned (owns Hamiltonian + accelerator)
    std::unique_ptr<qchem::ChargeDensity::DM_CD> itsDensity;
    std::vector<occ_t>                   itsOccupied;
    Observer                             itsObserver;
};

} // namespace qchem
