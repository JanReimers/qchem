// File: Calculation/Imp/Calculation.C  Implementation of the qchem::Calculation facade.
module;
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <set>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <nlohmann/json.hpp>
module qchem.Calculation;

import qchem.BasisSet.Molecule.Factory;             // BasisSet::Molecule::Factory
import qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet; // SymmetryAdaptedBasisSet (return of SymmetryAdapt)
import qchem.BasisSet.Molecule.PG_Cart.SymmetryAdapt;   // PG::SymmetryAdapt (the SALC builder)
import qchem.ElectronConfiguration.Molecule;        // Molecule_EC (global aufbau)
import qchem.PeriodicTable;                          // thePeriodicTable (Z -> element symbol, for the PP lookup)
import qchem.Pseudopotential.GTH_Potentials;         // GetGTH (Zion for the valence electron count + PP model)
import qchem.SCFAccelerator.Factory;                // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                           // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                               // Orbital, Orbitals

namespace qchem
{

namespace PG = ::qchem::BasisSet::Molecule::PG_Cart;
using SCFIter = qchem::SCFIterator::SCFIterator;

// Build the molecular orbital basis, optionally SALC point-group-blocked.  When symmetry is off the
// caller owns the raw basis directly; when on, the raw basis is handed to PG::SymmetryAdapt as a
// shared_ptr that the returned SymmetryAdaptedBasisSet keep-alives -- so the facade owns (and deletes)
// only the adapted basis, and the raw lifetime rides along.  PG::SymmetryAdapt throws if the basis has
// no Cartesian PGData IBS (the guard); the facade only ever builds that basis today, so it is safe.
static BasisSet::Real_BS* BuildBasis(const CalcOptions& opts, const std::shared_ptr<const Structure>& st)
{
    const bool spherical = (opts.angular == Angular::Spherical);
    // SALC now supports the Cartesian PGData basis AND the in-house MnD-spherical basis (SphData;
    // doc/SphericalSALCPlan.md S3a/S4).  Only libcint-spherical is unwired (S3b) -- and it presents as a
    // PGData with spherical components, which the Cartesian extractor would silently misread, so reject that
    // one combination clearly here rather than let SymmetryAdapt take the wrong path.
    if (opts.symmetry && spherical && opts.engine == Engine::LibCint)
        throw std::runtime_error("qchem::Calculation: {.symmetry=true} with angular=Spherical requires "
                                 "engine=MnD (libcint-spherical SALC is not yet supported)");

    const char* engine  = (opts.engine  == Engine::LibCint) ? "libcint" : "mnd";
    const char* angular = spherical ? "spherical" : "cartesian";
    BasisSet::Real_BS* raw = BasisSet::Molecule::Factory(
        nlohmann::json{{"basis", opts.basis}, {"engine", engine}, {"angular", angular}}, st.get());
    if (!opts.symmetry) return raw;
    std::shared_ptr<const BasisSet::Real_BS> rawShared(raw);
    return PG::SymmetryAdapt(rawShared, *st, opts.symmetryTol);
}

// Build the molecular electron configuration from the requested multiplicity 2S+1.  multiplicity<=0 is the
// minimal-spin default (closed-shell singlet for even Ne, doublet for odd, via Molecule_EC(Ne)); an explicit
// value is converted to (nUp,nDown) and parity-validated against Ne (fail loud, not a silent miscount).
static Molecule_EC* MakeMoleculeEC(int Ne, int multiplicity)
{
    if (multiplicity <= 0) return new Molecule_EC(Ne);    // auto / minimal spin (historical behaviour)
    const int twoS = multiplicity - 1;                    // = nUp - nDown
    if (twoS > Ne || ((Ne - twoS) % 2) != 0)
        throw std::runtime_error("qchem::Calculation: multiplicity " + std::to_string(multiplicity) +
            " (2S=" + std::to_string(twoS) + ") is incompatible with " + std::to_string(Ne) +
            " electrons -- need 0 <= 2S <= Ne and Ne-2S even.");
    return new Molecule_EC((Ne + twoS) / 2, (Ne - twoS) / 2);   // (nUp, nDown)
}

// The pseudopotential valence charge (Zion) for element \a Z, from its default-valence GTH entry.
static int PPZion(int Z) {return Pseudopotential::GetGTH(thePeriodicTable().GetSymbol(Z)).zion;}

// Re-express a structure as its VALENCE ions for a pseudopotential run: each atom keeps its true species Z
// (so the PP lookup and the true geometry are intact) but carries net charge Z-Zion(its own element), so
// GetNumElectrons() reports the total Zion valence count the EC + FittedVee charge constraint consume.
static std::shared_ptr<Structure> MakeValenceStructure(const Structure& st)
{
    auto mol = std::make_shared<Molecule>();
    for (size_t a=0; a<st.GetNumAtoms(); a++)
    {
        const int Z = st[a]->itsZ;
        mol->Insert(new Atom(Z, double(Z - PPZion(Z)), st[a]->itsR));   // per-element charge = Z - Zion
    }
    return mol;
}

// The distinct pseudopotential species (element, valence) present in \a st -- the per-Z router the
// multi-species PP Hamiltonian is built from.  A single-species molecule yields a 1-element list.
static std::vector<std::pair<std::string,int>> PPSpecies(const Structure& st)
{
    std::vector<std::pair<std::string,int>> species;
    std::set<int> seen;
    for (size_t a=0; a<st.GetNumAtoms(); a++)
    {
        const int Z = st[a]->itsZ;
        if (seen.insert(Z).second) species.push_back({thePeriodicTable().GetSymbol(Z), PPZion(Z)});
    }
    return species;
}

Calculation::Calculation(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc)
    : itsStructure(std::make_shared<Molecule>(st))   // deep copy: the facade owns its own structure
    , itsOpts(opts)
    , itsAcc(acc)
{
    // An open shell (2S>0) is unrestricted: it needs distinct up/down densities, so promote to polarized.
    if (itsOpts.multiplicity > 1) itsOpts.pol = Pol::Polarized;

    // A pseudopotential run works on the VALENCE ions (Z-Zion charge) so the electron count is the valence
    // count -- built before the basis/EC, which both read it off the structure.
    if (itsOpts.pseudopotential) itsStructure = MakeValenceStructure(*itsStructure);

    itsBasis = BuildBasis(itsOpts, itsStructure);                  // raw, or SALC-blocked if .symmetry
    itsEC    = MakeMoleculeEC(int(itsStructure->GetNumElectrons()), itsOpts.multiplicity);
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
    // orbital basis, xalpha) are ignored for HF/1-e/Dirac.  A pseudopotential run takes the PP front door
    // instead (LSDA valence Hamiltonian: V_loc + KB projectors + Zion ion-ion in place of Ven).
    auto* ham = itsOpts.pseudopotential
        ? qchem::Hamiltonian::Factory(itsStructure, PPSpecies(*itsStructure), itsOpts.mesh, itsBasis)
        : qchem::Hamiltonian::Factory(itsOpts.model, itsOpts.pol, itsStructure,
                                      itsOpts.mesh, itsBasis, itsOpts.xalpha);

    // DFT (and the LSDA pseudopotential) need DIIS engaged from iteration 0: the default Z-scaled EMax gate
    // keeps DIIS off and the SCF limit-cycles to a non-converged snapshot (the "M_Sym layout UB" mechanism,
    // see M_DFT).  So default them to a high EMax ("DIIS from start") unless the caller pinned one.
    const bool dftLike = dft || itsOpts.pseudopotential;
    const double emax = itsAcc.eMax > 0.0 ? itsAcc.eMax : (dftLike ? 100.0 : Z*Z*0.1/32);
    nlohmann::json jsacc = {{"NProj", itsAcc.nProj}, {"EMax", emax},
                            {"EMin", itsAcc.eMin}, {"SVTol", itsAcc.svTol}};
    auto* accel = qchem::SCFAccelerators::Factory(qchem::SCFAccelerators::Type::DIIS, jsacc);

    delete itsScf;                                     // releases the previous Hamiltonian + accelerator
    // Seed: an explicit opts.seed wins; otherwise auto -- DFT from superposition-of-atomic-densities
    // (SAD), HF/1-e from the core guess (Default).
    using qchem::ChargeDensity::SeedStrategy;
    const auto seed = (itsOpts.seed != SeedStrategy::Default) ? itsOpts.seed
                      : (dft ? SeedStrategy::SAD : SeedStrategy::Default);
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
