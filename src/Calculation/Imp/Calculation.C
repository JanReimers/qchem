// File: Calculation/Imp/Calculation.C  Implementation of the qchem::Calculation facade.
module;
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cassert>
#include <nlohmann/json.hpp>
module qchem.Calculation;

import qchem.BasisSet.Molecule.Factory;      // BasisSet::Molecule::Factory
import qchem.ElectronConfiguration.Molecule; // Molecule_EC (global aufbau)
import qchem.SCFAccelerator.Factory;         // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                    // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                        // Orbital, Orbitals

namespace qchem
{

using SCFIter = qchem::SCFIterator::SCFIterator;

Calculation::Calculation(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc)
    : itsStructure(std::make_shared<Molecule>(st))   // deep copy: the facade owns its own structure
    , itsOpts(opts)
    , itsAcc(acc)
{
    itsBasis = BasisSet::Molecule::Factory(nlohmann::json{{"basis", opts.basis}}, itsStructure.get());
    itsEC    = new Molecule_EC(int(itsStructure->GetNumElectrons()));   // aufbau, global across irreps
    Converge();                                       // so Energy()/Density() are ready on return
}

Calculation::~Calculation()
{
    delete itsScf;     // owns + deletes the Hamiltonian and accelerator
    delete itsBasis;
    delete itsEC;
}

bool Calculation::Converge(const SCFParams& params)
{
    const int Z = itsStructure->GetNuclearCharge();
    auto* ham = qchem::Hamiltonian::Factory(itsOpts.model, itsOpts.pol, itsStructure);

    nlohmann::json jsacc = {{"NProj", itsAcc.nProj},
                            {"EMax",  itsAcc.eMax > 0.0 ? itsAcc.eMax : Z*Z*0.1/32},
                            {"EMin",  itsAcc.eMin},
                            {"SVTol", itsAcc.svTol}};
    auto* accel = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);

    delete itsScf;                                     // releases the previous Hamiltonian + accelerator
    itsScf = new SCFIter(itsBasis, itsEC, ham, accel,
                         qchem::ChargeDensity::SeedStrategy::Default, itsStructure.get());
    if (itsObserver) itsScf->SetObserver(itsObserver);

    bool ok = itsScf->Iterate(params);
    RebuildSampling();
    return ok;
}

void Calculation::RebuildSampling()
{
    const auto* wf = itsScf->GetWaveFunction();
    itsDensity.reset(wf->GetChargeDensity());          // GetChargeDensity hands back a heap rho -- we own it

    itsOccupied.clear();
    for (const Irrep& irr : wf->GetQNs())
        for (const auto* o : wf->GetOrbitals(irr)->Iterate())   // Iterate() yields const Orbital*
            if (o->IsOccupied())
            {
                // An Orbital IS a ScalarFunction<double> (cast between abstract faces, per CLAUDE.md).
                const auto* phi = dynamic_cast<const sf_t*>(o);
                assert(phi);
                itsOccupied.push_back({o->GetEigenEnergy(), phi});
            }
    std::sort(itsOccupied.begin(), itsOccupied.end(),
              [](const occ_t& a, const occ_t& b){ return a.first < b.first; });   // ascending energy
}

void Calculation::OnIteration(Observer obs) {itsObserver = std::move(obs);}

double                 Calculation::Energy()      const {return itsScf->GetEnergy().GetTotalEnergy();}
qchem::EnergyBreakdown Calculation::EnergyTerms() const {return itsScf->GetEnergy();}

const Calculation::sf_t& Calculation::Density() const
{
    assert(itsDensity);
    return *itsDensity;
}

const Calculation::sf_t& Calculation::HOMO() const
{
    assert(!itsOccupied.empty());
    return *itsOccupied.back().second;
}

const Calculation::sf_t& Calculation::Orbital(size_t i) const
{
    assert(i < itsOccupied.size());
    return *itsOccupied[i].second;
}

size_t Calculation::NumOccupied() const {return itsOccupied.size();}

const Calculation::orbitals_t* Calculation::Orbitals(const Irrep& qns) const
{
    return itsScf->GetWaveFunction()->GetOrbitals(qns);
}

size_t Calculation::IterationCount() const {return itsScf->GetIterationCount();}
bool   Calculation::IsConverged()    const {return itsScf->Converged();}

} // namespace qchem
