// File: Calculation/Imp/SolidCalculation.C  The periodic SCF front door.
//
// Imports ZERO .Internal. modules -- the same discipline Imp/Calculation.C and Imp/AtomCalculation.C
// already keep, and the reason Step 4 had to open two public doors before this file could exist.
module;
#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
module qchem.SolidCalculation;

import qchem.BasisSet.Molecule.LatticeSum1E;   // LatticeSum1E::MaxExponent -- alpha_max (a capability face)
import qchem.BasisSet.Orbital_1E_IBS;          // Real_OIBS / Complex_OIBS (the per-irrep bases to iterate)
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH -> HGH local PP -> alpha_pp
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.Mesh.XCPolicy;                    // XCMeshSharpness / ResolveXCMesh (the grid decision)
import qchem.ElectronConfiguration.Crystal;    // Crystal_EC
import qchem.Symmetry.Irrep;                   // Irrep, Spin

namespace qchem
{

//---------------------------------------------------------------------------------------------------
//  The run's SHARPNESS -- an above-SCFIterator decision input.
//
//  Both sources come off ABSTRACT capability faces via the sanctioned abstract->abstract cross-cast, so
//  nothing here touches a concrete basis or a concrete PP model:
//    alpha_max -- BasisSet::Molecule::LatticeSum1E::MaxExponent(), documented there as "the GPW
//                 density-grid cutoff floor".
//    alpha_pp  -- Pseudopotential::LocalPotential_Gaussian::ShortRangeGaussian(Z), whose terms carry
//                 alpha = 1/(2 r_loc^2).  A model with no closed-Gaussian short part does not implement
//                 the face; that leaves alpha_pp at 0, which the selector reads as "not measurable" --
//                 NOT as "smooth".
//---------------------------------------------------------------------------------------------------
static qcMesh::XCMeshSharpness GatherSharpness(const Lattice_3D& lat, const BasisSet::Real_BS& mol,
                                               const SolidCalcOptions& o)
{
    qcMesh::XCMeshSharpness s;
    s.cellEdge = lat.GetUnitCell().GetMaximumCellEdge();
    s.nAtoms   = int(lat.GetUnitCell().GetNumAtoms());
    s.imposed  = o.imposeSymmetry;
    for (auto ibs : const_cast<BasisSet::Real_BS&>(mol).Iterate<BasisSet::Real_OIBS>())
        if (const auto* ls=dynamic_cast<const BasisSet::Molecule::LatticeSum1E*>(ibs))
            { s.alphaMax = ls->MaxExponent(); break; }
    for (const auto& [element, valence] : o.species)
    {
        const int Z = int(thePeriodicTable().GetZ(element));
        const Pseudopotential::HGH_LocalPotential loc = Pseudopotential::GetGTH(element,"LDA",valence).local;
        const auto* g = static_cast<const Pseudopotential::LocalPotential_Gaussian*>(&loc);
        for (const auto& t : g->ShortRangeGaussian(Z)) s.alphaPP = std::max(s.alphaPP, t.alpha);
    }
    return s;
}

//---------------------------------------------------------------------------------------------------
struct SolidCalculation::Imp
{
    SolidCalcOptions opts;
    std::shared_ptr<const Structure>            st;
    std::unique_ptr<BasisSet::Complex_BS>       bs;
    std::unique_ptr<Crystal_EC>                 ec;
    qchem::Hamiltonian::cHamiltonian*           ham = nullptr;   // owned by the iterator once handed over
    qcMesh::MeshParams                          xcMesh;          // AFTER Auto resolution
    std::unique_ptr<qchem::SCFIterator::SolidSCFIterator> scf;
    std::unique_ptr<qchem::ChargeDensity::cDM_CD>         cd;    // the converged density (outlives the WF)
    SCFAccelerators::SolidAcceleratorOptions    accOpts;
    bool   converged = false;
    double charge    = 0.0;
};

//---------------------------------------------------------------------------------------------------
SolidCalculation::SolidCalculation(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                   const SolidCalcOptions& opts, const SCFParams& params,
                                   const SCFAccelerators::SolidAcceleratorOptions& acc)
    : itsImp(std::make_unique<Imp>())
{
    itsImp->opts    = opts;
    itsImp->accOpts = acc;
    itsImp->st      = lat.GetStructure();

    namespace L3 = BasisSet::Lattice_3D;
    // THE WORKING-TYPE DECISION (doc/RealComplexPlan.md §3, Step 3c-3): a block is real ⇔ its irrep is
    // (TRIM) ∧ every term preserves realness.  This composition root builds the LDA GPW stack --
    // kinetic, the PP trio, Hartree, XC -- every member of which PreservesReal(), so the term half is
    // TRUE here; it is asserted against the constructed Hamiltonian below, so a future term that
    // breaks realness must also flip this forecast.  forceComplex is the §6 ansatz-policy downgrade.
    const bool hamPreservesReal = !opts.forceComplex;
    itsImp->bs.reset(L3::GPWFactory(lat, mol, L3::GPWParams{
        .densityEcut = opts.densityEcut, .cutoffFactor = opts.cutoffFactor, .raster = opts.raster,
        .images = opts.images, .kShift = opts.kShift, .ladderFactor = opts.ladderFactor,
        .imposeSymmetry = opts.imposeSymmetry, .hamPreservesReal = hamPreservesReal}));

    // DECISION 1 -- the XC quadrature.  Resolve Auto HERE, once, from facts about the run.  Downstream
    // consumers compare ==Becke, so an unresolved Auto would silently read as Uniform; resolving it at the
    // point the spec enters the Hamiltonian is what makes that impossible.
    itsImp->xcMesh = qcMesh::ResolveXCMesh(opts.xcMesh, GatherSharpness(lat, *mol, opts));

    // DECISION 2 -- the spin bookkeeping.  multiplicity -> (nUp,nDown), with the parity check that catches
    // a singlet asked of an odd electron count BEFORE integer division silently empties a channel.
    const int twoS = opts.multiplicity>1 ? opts.multiplicity-1 : opts.Nelec%2;
    const bool polarized = opts.multiplicity>=1;
    if ((opts.Nelec-twoS)%2!=0 || twoS>opts.Nelec)
        throw std::runtime_error("SolidCalcOptions: multiplicity "+std::to_string(opts.multiplicity)
                                 +" parity disagrees with Nelec "+std::to_string(opts.Nelec));
    auto irreps = itsImp->bs->GetIrreps(Spin::None);   // one Bloch irrep per BZ k-block (weights carry Sum_k)
    itsImp->ec = std::make_unique<Crystal_EC>(irreps, (opts.Nelec+twoS)/2, (opts.Nelec-twoS)/2,
                                              opts.globalFermi);

    itsImp->ham = qchem::Hamiltonian::Factory(
        polarized ? qchem::Hamiltonian::Pol::Polarized : qchem::Hamiltonian::Pol::UnPolarized,
        itsImp->st, itsImp->bs.get(), opts.species, "LDA", itsImp->xcMesh, opts.vxcFit);
    // The forecast crosscheck: the basis was built on the promise that every term preserves realness
    // (the AND's term half, above); the constructed Hamiltonian must agree, or real blocks were built
    // that its terms cannot serve.
    assert((!hamPreservesReal || itsImp->ham->PreservesReal()) &&
           "SolidCalculation: the term stack no longer preserves realness -- update the forecast above");

    // DECISION 3 -- the accelerator, by policy, through the public typed door.
    auto* accel = SCFAccelerators::Factory(opts.accelerator, acc);

    itsImp->scf = std::make_unique<qchem::SCFIterator::SolidSCFIterator>(
        itsImp->bs.get(), itsImp->ec.get(), itsImp->ham, accel,
        opts.seed, itsImp->st.get(), opts.ortho, opts.orthoTol);

    (void)Converge(params);   // the ctor ATTEMPTS; the caller faces the result via Result()
}

SolidCalculation::~SolidCalculation() = default;

//---------------------------------------------------------------------------------------------------
// THE ANSWERS LIVE ON THE PROOF (doc/OpenWork.md N1/T1).  Each of these used to sit on SolidCalculation
// itself with no precondition, so a run that never converged served a plausible number to anyone who did
// not think to ask Converged() first.  Reaching them now REQUIRES having been handed a Converged, and the
// only source of one is a successful attempt.
double SolidCalculation::Converged::Energy()      const {return itsImp->scf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown SolidCalculation::Converged::EnergyTerms() const {return itsImp->scf->GetEnergy();}
double SolidCalculation::Converged::TotalCharge() const {return itsImp->charge;}
size_t SolidCalculation::Converged::IterationCount() const {return itsImp->scf->GetIterationCount();}

const ScalarFunction<double>& SolidCalculation::Converged::Density() const
{
    // Still an assert, and legitimately so: reaching here without a density would be a BROKEN INVARIANT
    // (a Converged is only ever minted beside one), not a user error.  The user-error case is what the
    // type now prevents outright.
    assert(itsImp->cd && "SolidCalculation::Converged::Density: converged without a density");
    return *itsImp->cd;
}

//---------------------------------------------------------------------------------------------------
//! Mint the outcome from the state the last attempt left behind -- one place, so Converge() and Result()
//! cannot drift apart in what they call a success.
Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Outcome_() const
{
    using O = Outcome<Converged, SCFFailure>;
    if (itsImp->converged) return O::Ok(Converged(itsImp.get()));
    SCFFailure f;
    f.why        = SCFFailure::Why::NotConverged;
    f.iterations = itsImp->scf->GetIterationCount();
    f.lastEnergy = itsImp->scf->GetEnergy().GetTotalEnergy();
    f.details    = "the SCF reached its iteration limit with the residual still above tolerance";
    return O::Fail(std::move(f));
}

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Converge(const SCFParams& params)
{
    assert(itsImp->scf);
    itsImp->scf->Iterate(params);
    itsImp->converged = itsImp->scf->Converged();
    // Take the density OUT of the wave function: it stays valid after the iterator's state moves on
    // because its basis block is `bs`, which this object owns and outlives it.
    auto cd = itsImp->scf->GetWaveFunction()->GetChargeDensity();   // BUILT for us; we take it
    itsImp->charge = cd->GetTotalCharge();
    itsImp->cd = std::move(cd);
    return Outcome_();
}

Outcome<SolidCalculation::Converged, SCFFailure> SolidCalculation::Result() const {return Outcome_();}

qchem::EnergyBreakdown SolidCalculation::LastIterateTerms()  const {return itsImp->scf->GetEnergy();}
double                 SolidCalculation::LastIterateCharge() const {return itsImp->charge;}

void   SolidCalculation::OnIteration(Observer obs)      { itsImp->scf->SetObserver(std::move(obs)); }
bool   SolidCalculation::DidConverge()    const         { return itsImp->converged; }
size_t SolidCalculation::IterationCount() const         { return itsImp->scf->GetIterationCount(); }

const qcMesh::MeshParams& SolidCalculation::ResolvedXCMesh() const { return itsImp->xcMesh; }

} //namespace qchem
