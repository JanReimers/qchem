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

import qchem.BasisSet.Molecule.Factory;             // BasisSet::Molecule::Factory
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet; // SymmetryAdaptedBasisSet (return of SymmetryAdapt)
import qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt;   // PG::SymmetryAdapt (the SALC builder)
import qchem.ElectronConfiguration.Molecule;        // Molecule_EC (global aufbau)
import qchem.SCFAccelerator.Factory;                // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                           // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                               // Orbital, Orbitals

namespace qchem
{

namespace PG = ::BasisSet::Molecule::PG_Cart;
using SCFIter = qchem::SCFIterator::SCFIterator;

// Build the molecular orbital basis, optionally SALC point-group-blocked.  When symmetry is off the
// caller owns the raw basis directly; when on, the raw basis is handed to PG::SymmetryAdapt as a
// shared_ptr that the returned SymmetryAdaptedBasisSet keep-alives -- so the facade owns (and deletes)
// only the adapted basis, and the raw lifetime rides along.  PG::SymmetryAdapt throws if the basis has
// no Cartesian PGData IBS (the guard); the facade only ever builds that basis today, so it is safe.
static BasisSet::Real_BS* BuildBasis(const CalcOptions& opts, const std::shared_ptr<const Structure>& st)
{
    BasisSet::Real_BS* raw = BasisSet::Molecule::Factory(nlohmann::json{{"basis", opts.basis}}, st.get());
    if (!opts.symmetry) return raw;
    std::shared_ptr<const BasisSet::Real_BS> rawShared(raw);
    return PG::SymmetryAdapt(rawShared, *st, opts.symmetryTol);
}

Calculation::Calculation(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc)
    : itsStructure(std::make_shared<Molecule>(st))   // deep copy: the facade owns its own structure
    , itsOpts(opts)
    , itsAcc(acc)
{
    itsBasis = BuildBasis(opts, itsStructure);                         // raw, or SALC-blocked if .symmetry
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
    const int  Z   = itsStructure->GetNuclearCharge();
    const bool dft = qchem::Hamiltonian::IsDFT(itsOpts.model);

    // The unified resolver turns the Model token into the concrete Hamiltonian; the DFT extras (mesh,
    // orbital basis, xalpha) are ignored for HF/1-e/Dirac.
    auto* ham = qchem::Hamiltonian::Factory(itsOpts.model, itsOpts.pol, itsStructure,
                                            itsOpts.mesh, itsBasis, itsOpts.xalpha);

    // DFT needs DIIS engaged from iteration 0: the default Z-scaled EMax gate keeps DIIS off and the DFT
    // SCF limit-cycles to a non-converged snapshot (the "M_Sym layout UB" mechanism, see M_DFT).  So
    // default DFT to a high EMax ("DIIS from start") unless the caller pinned one explicitly.
    const double emax = itsAcc.eMax > 0.0 ? itsAcc.eMax : (dft ? 100.0 : Z*Z*0.1/32);
    nlohmann::json jsacc = {{"NProj", itsAcc.nProj}, {"EMax", emax},
                            {"EMin", itsAcc.eMin}, {"SVTol", itsAcc.svTol}};
    auto* accel = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);

    delete itsScf;                                     // releases the previous Hamiltonian + accelerator
    // DFT seeds from superposition-of-atomic-densities (SAD); HF/1-e take the core guess (Default).
    const auto seed = dft ? qchem::ChargeDensity::SeedStrategy::SAD
                          : qchem::ChargeDensity::SeedStrategy::Default;
    itsScf = new SCFIter(itsBasis, itsEC, ham, accel, seed, itsStructure.get());
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
