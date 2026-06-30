// File: Calculation/Imp/AtomCalculation.C  Implementation of the qchem::AtomCalculation facade.
module;
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cassert>
#include <nlohmann/json.hpp>
module qchem.AtomCalculation;

import qchem.Hamiltonian.Factory;                   // Factory(Model,Pol,st,mesh,bs,xalpha), IsDFT
import qchem.ElectronConfiguration.AtomNR;          // Atom_EC
import qchem.ElectronConfiguration.AtomDirac;       // AtomDirac_EC
import qchem.SCFAccelerator.Factory;                // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                           // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                               // Orbital, Orbitals
import qchem.Types;                                  // Vector3D / rvec3_t

namespace qchem
{

using SCFIter = qchem::SCFIterator::SCFIterator;

static bool IsDirac(Model m) {return m == Model::DE1 || m == Model::DHF;}

// Build the atomic exponent-pool basis.  An explicit pool ({type,N,emin,emax}, opts.N>0) reproduces the
// A_SG/A_SL exponent-sweep json the scaffold built by hand; otherwise the BasisSetAccuracy preset is used
// (the EC fixes which angular irreps the pool spans).
static BasisSet::Real_BS* BuildBasis(const AtomCalcOptions& opts, int Z, const ElectronConfiguration& ec)
{
    if (opts.N > 0)
        return BasisSet::Atom::Factory(nlohmann::json{{"type", opts.type}, {"N", opts.N},
                                                      {"emin", opts.emin}, {"emax", opts.emax}}, Z);
    return BasisSet::Atom::Factory(opts.accuracy, opts.type, Z, ec);
}

AtomCalculation::AtomCalculation(int Z, int charge, const AtomCalcOptions& opts, const AcceleratorOptions& acc)
    : itsZ(Z)
    , itsNe(Z - charge)
    , itsOpts(opts)
    , itsAcc(acc)
    , itsStructure(std::make_shared<Atom>(Z, double(charge), rvec3_t(0,0,0)))
{
    // The EC follows the model: relativistic models need the Dirac (j-coupled) configuration.
    itsEC = IsDirac(opts.model) ? static_cast<ElectronConfiguration*>(new AtomDirac_EC(itsNe))
                                : static_cast<ElectronConfiguration*>(new Atom_EC(itsNe));
    itsBasis = BuildBasis(opts, Z, *itsEC);
    Converge();
}

AtomCalculation::~AtomCalculation()
{
    delete itsScf;     // owns + deletes the Hamiltonian and accelerator
    delete itsBasis;
    delete itsEC;
}

bool AtomCalculation::Converge(const SCFParams& params)
{
    auto* ham = qchem::Hamiltonian::Factory(itsOpts.model, itsOpts.pol, itsStructure,
                                            itsOpts.mesh, itsBasis, itsOpts.xalpha);

    // Atoms use the proven Z-scaled DIIS gate (EMax = Z^2*0.1/32) unless the caller pins one.
    const double emax = itsAcc.eMax > 0.0 ? itsAcc.eMax : itsZ*itsZ*0.1/32;
    nlohmann::json jsacc = {{"NProj", itsAcc.nProj}, {"EMax", emax},
                            {"EMin", itsAcc.eMin}, {"SVTol", itsAcc.svTol}};
    auto* accel = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);

    delete itsScf;
    // Default seed for atoms is the core guess (atoms never use the molecular SAD seed).
    using qchem::ChargeDensity::SeedStrategy;
    const auto seed = (itsOpts.seed != SeedStrategy::Default) ? itsOpts.seed : SeedStrategy::CoreGuess;
    itsScf = new SCFIter(itsBasis, itsEC, ham, accel, seed, itsStructure.get());
    if (itsObserver) itsScf->SetObserver(itsObserver);

    bool ok = itsScf->Iterate(params);
    RebuildSampling();
    return ok;
}

void AtomCalculation::RebuildSampling()
{
    const auto* wf = itsScf->GetWaveFunction();
    itsDensity.reset(wf->GetChargeDensity());

    itsOccupied.clear();
    for (const Irrep& irr : wf->GetQNs())
        for (const auto* o : wf->GetOrbitals(irr)->Iterate())
            if (o->IsOccupied())
            {
                const auto* phi = dynamic_cast<const sf_t*>(o);
                assert(phi);
                itsOccupied.push_back({o->GetEigenEnergy(), phi});
            }
    std::sort(itsOccupied.begin(), itsOccupied.end(),
              [](const occ_t& a, const occ_t& b){ return a.first < b.first; });
}

void AtomCalculation::OnIteration(Observer obs) {itsObserver = std::move(obs);}

double                 AtomCalculation::Energy()      const {return itsScf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown AtomCalculation::EnergyTerms() const {return itsScf->GetEnergy();}
double                 AtomCalculation::TotalCharge() const
{
    return itsScf->GetWaveFunction()->GetChargeDensity()->GetTotalCharge();
}

const AtomCalculation::sf_t& AtomCalculation::Density() const {assert(itsDensity); return *itsDensity;}
const AtomCalculation::sf_t& AtomCalculation::HOMO()    const {assert(!itsOccupied.empty()); return *itsOccupied.back().second;}
const AtomCalculation::sf_t& AtomCalculation::Orbital(size_t i) const {assert(i<itsOccupied.size()); return *itsOccupied[i].second;}
size_t                       AtomCalculation::NumOccupied() const {return itsOccupied.size();}

const AtomCalculation::orbitals_t* AtomCalculation::Orbitals(const Irrep& qns) const
{
    return itsScf->GetWaveFunction()->GetOrbitals(qns);
}

size_t AtomCalculation::IterationCount() const {return itsScf->GetIterationCount();}
bool   AtomCalculation::IsConverged()    const {return itsScf->Converged();}

} // namespace qchem
