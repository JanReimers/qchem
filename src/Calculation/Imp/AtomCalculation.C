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

import qchem.Hamiltonian.Factory;                   // Factory(...), IsDFT, XCFunctional (XC selector)
import qchem.ElectronConfiguration.AtomNR;          // Atom_EC, PseudoAtom_EC
import qchem.ElectronConfiguration.AtomDirac;       // AtomDirac_EC
import qchem.PeriodicTable;                          // thePeriodicTable().GetSymbol(Z) (the PP element)
import qchem.SCFAccelerator.Factory;                // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                           // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                               // Orbital, Orbitals
import qchem.Types;                                  // Vector3D / rvec3_t
import qchem.BasisSet.IrrepBasisSet;                 // theCache<T>() / EmitReport (the run's `cache` section)
import qchem.Reporting;                              // report:: run sink + the scf-section writer

namespace qchem
{

using SCFIter = qchem::SCFIterator::MolecularSCFIterator;   // item 2: the atom/molecule display subclass

static bool IsDirac(Model m) {return m == Model::DE1 || m == Model::DHF;}

// Human-readable model tag for the run report's `scf` section (a PP run reports "PP" -- it bypasses the Model enum).
static std::string ModelName(Model m)
{
    switch (m)
    {
        case Model::HF:     return "HF";     case Model::E1:  return "E1";
        case Model::DHF:    return "DHF";    case Model::DE1: return "DE1";
        case Model::Xalpha: return "Xalpha"; case Model::LDA: return "LDA";
    }
    return "?";
}

// Atomic-basis family tag for the report (mirrors the scfrun --basis names).
static std::string AtomTypeName(AtomType t)
{
    switch (t)
    {
        case AtomType::Slater:       return "Slater";       case AtomType::Gaussian:     return "Gaussian";
        case AtomType::BSpline6:     return "BSpline6";      case AtomType::BSpliner6:    return "BSpliner6";
        case AtomType::Slater_RKB:   return "Slater_RKB";    case AtomType::Gaussian_RKB: return "Gaussian_RKB";
    }
    return "?";
}

// Build the run report's `scf` section (the atom-side analogue of Calculation.C's EmitScfSection): the
// "what recipe am I running" header.  The report IS json, so this provider owns the keys; standard = the
// settled recipe surface (doc/RunReportPlan.md schema).  A no-op console-wise unless a client attached one.
static void EmitScfSection(const AtomCalcOptions& o, int ne, const std::string& accelerator,
                           const SCFParams& p)
{
    namespace rpt = qchem::report;
    rpt::json standard;
    standard["model"]       = o.pseudopotential ? std::string("PP") : ModelName(o.model);
    standard["pol"]         = (o.pol == Pol::Polarized) ? "P" : "U";
    standard["nElectrons"]  = ne;
    standard["nMaxIter"]    = long(p.NMaxIter);
    standard["minDrho"]     = p.MinΔρ;
    standard["minFD"]       = p.MinFD;
    standard["minDE"]       = p.MinΔE;
    standard["accelerator"] = accelerator;
    if (o.pseudopotential) standard["valence"] = o.valence;

    rpt::json scf;
    scf["standard"] = standard;
    rpt::EmitSection("scf", scf);
}

// PP Zion: the explicit valence (Zion) if given, else the electron count (a neutral pseudo-atom).
static int PPZion(const AtomCalcOptions& opts, int Ne) {return opts.valence > 0 ? opts.valence : Ne;}

// Build the atomic exponent-pool basis.  An explicit pool ({type,N,emin,emax}, opts.N>0) reproduces the
// A_SG/A_SL exponent-sweep json the scaffold built by hand; otherwise the BasisSetAccuracy preset is used
// (the EC fixes which angular irreps the pool spans).
static BasisSet::Real_BS* BuildBasis(const AtomCalcOptions& opts, int Z, const ElectronConfiguration& ec)
{
    // Explicit pool: build from the EC (which fixes the angular irreps), NOT the nuclear Z -- the
    // Factory(json,Z) overload treats Z as an electron count, so for an ion that would over-populate the
    // angular irreps (e.g. a 1-electron hydrogenic ion would gain p/d shells).
    // Independent per-l exponent lists -- highest precedence (the generator's per-l validation path).
    if (!opts.exponentsByL.empty())
    {
        nlohmann::json shells=nlohmann::json::array();
        for (const auto& [l,es] : opts.exponentsByL) shells.push_back({{"l",l},{"exponents",es}});
        return BasisSet::Atom::Factory(nlohmann::json{{"type", opts.type}, {"shells", shells}}, ec);
    }
    // Explicit single exponent list ("bring your own exponents") -- applied to every occupied l (with ltrim).
    if (!opts.exponents.empty())
        return BasisSet::Atom::Factory(nlohmann::json{{"type", opts.type},
                                                      {"exponents", opts.exponents}, {"ltrim", opts.ltrim}}, ec);
    if (opts.N > 0)
        return BasisSet::Atom::Factory(nlohmann::json{{"type", opts.type}, {"N", opts.N},
                                                      {"emin", opts.emin}, {"emax", opts.emax}}, ec);
    return BasisSet::Atom::Factory(opts.accuracy, opts.type, Z, ec);
}

AtomCalculation::AtomCalculation(int Z, int charge, const AtomCalcOptions& opts,
                                 const SCFParams& params, const AcceleratorOptions& acc)
    : itsZ(Z)
    , itsNe(Z - charge)
    , itsOpts(opts)
    , itsAcc(acc)
    , itsStructure(std::make_shared<Atom>(Z, double(charge), rvec3_t(0,0,0)))
{
    // The EC follows the model: a pseudo-atom uses the valence-only PseudoAtom_EC (keyed by the full Z and
    // the net charge Zion-electrons); relativistic models need the Dirac (j-coupled) configuration;
    // everything else is a plain Atom_EC.
    itsEC = opts.pseudopotential ? static_cast<ElectronConfiguration*>(new PseudoAtom_EC(itsZ, PPZion(opts,itsNe)-itsNe))
          : IsDirac(opts.model)  ? static_cast<ElectronConfiguration*>(new AtomDirac_EC(itsNe))
          :                        static_cast<ElectronConfiguration*>(new Atom_EC(itsNe));
    itsBasis = BuildBasis(opts, Z, *itsEC);
    Converge(params);
}

AtomCalculation::~AtomCalculation()
{
    delete itsScf;     // owns + deletes the Hamiltonian and accelerator
    delete itsBasis;
    delete itsEC;
}

bool AtomCalculation::Converge(const SCFParams& params)
{
    namespace H = qchem::Hamiltonian;
    // Three DFT/HF routes: a pseudopotential (PP front door, valence electrons = itsNe), an explicit XC
    // functional override (the public selector, e.g. libxc), or the model's built-in Hamiltonian/functional.
    auto* ham = itsOpts.pseudopotential
        ? H::Factory(itsOpts.pol, itsStructure, thePeriodicTable().GetSymbol(itsZ), PPZion(itsOpts,itsNe), itsOpts.mesh, itsBasis)
        : itsOpts.xc.has_value()
            ? H::Factory(itsOpts.pol, itsStructure, *itsOpts.xc, itsOpts.mesh, itsBasis)
            : H::Factory(itsOpts.model, itsOpts.pol, itsStructure, itsOpts.mesh, itsBasis, itsOpts.xalpha);

    // Atoms use the proven Z-scaled DIIS gate (EMax = Z^2*0.1/32) unless the caller pins one.  The json
    // escape hatch (opts.accelerator) then merges over these defaults and selects the accelerator TYPE --
    // the accelerator-tuning surface the scfrun driver needs (DIIS/GDM/Ladder/directmin).
    const double emax = itsAcc.eMax > 0.0 ? itsAcc.eMax : itsZ*itsZ*0.1/32;
    nlohmann::json jsacc = {{"NProj", itsAcc.nProj}, {"EMax", emax},
                            {"EMin", itsAcc.eMin}, {"SVTol", itsAcc.svTol}};
    std::string ts = "DIIS";
    if (itsOpts.accelerator.is_object() && !itsOpts.accelerator.empty())
    {
        for (auto& [k,v] : itsOpts.accelerator.items()) jsacc[k]=v;
        ts = itsOpts.accelerator.value("type", std::string("DIIS"));
    }
    using AType = qchem::SCFAccelerators::Type;
    const bool directmin = (ts=="directmin" || ts=="DirectMin");
    if (directmin) jsacc["EMax"]=1e10;   // GDM always steps; the direct-min loop seeds via diagonalize
    const AType atype = directmin ? AType::GDM : ts=="Ladder" ? AType::Ladder : ts=="GDM" ? AType::GDM : AType::DIIS;
    auto* accel = qchem::SCFAccelerators::Factory(atype, jsacc);

    delete itsScf;
    // Default seed for atoms is the core guess (atoms never use the molecular SAD seed).
    using qchem::ChargeDensity::SeedStrategy;
    const auto seed = (itsOpts.seed != SeedStrategy::Default) ? itsOpts.seed : SeedStrategy::CoreGuess;
    // RELATIVISTIC (Dirac) bases are even-tempered / kinetic-balanced: their overlap is near-singular BY
    // CONSTRUCTION, and canonical-ortho mode-dropping breaks kinetic balance -> variational collapse.  So a
    // Dirac SCF ALWAYS uses PLAIN Cholesky (never truncate); a non-relativistic atom uses opts.ortho (Auto by
    // default, but Cholesky when a caller deliberately wants an over-complete basis kept whole).
    const qchem::Ortho ortho = IsDirac(itsOpts.model) ? qchem::Cholesky : itsOpts.ortho;

    // Open a run report for the whole atomic SCF -- the atom-side mirror of qchem::Calculation (Reporting is
    // new; this brings the atomic front door onto it, so scfrun/valgen and the A_* tests all record a run).
    // The RAII guard guarantees End() even if a step throws, so a failed run never leaks the sink's context/
    // cursor stacks.  Nothing renders unless a client attached a console (report::SetConsole); nested runs (a
    // valence-basis generator brackets its own) stay console-quiet below depth 1.
    namespace rpt = qchem::report;
    rpt::Begin(itsStructure->ID());
    struct RunScope { ~RunScope() { rpt::End(); } } runScope;
    EmitScfSection(itsOpts, itsNe, ts, params);        // the "what recipe am I running" header

    // The `basis` section is ASSEMBLED across layers: stamp the scalars here, then the WF build (inside
    // `new SCFIter`) fills basis.perIrrep/removed via the cursor (the shared LASolver conditioning), and the
    // Section renders it once, complete, when this scope closes.
    {
        rpt::Section basis("basis");
        rpt::Set("type",       AtomTypeName(itsOpts.type));
        rpt::Set("nFunctions", (long)itsBasis->GetNumFunctions());
        itsScf = new SCFIter(itsBasis, itsEC, ham, accel, seed, itsStructure.get(), ortho);
    }
    // (directmin => AType::GDM above, whose WantsLineSearch() drives the direct-min loop -- no SetDirectMin needed.)
    if (itsObserver) itsScf->SetObserver(itsObserver);

    bool ok = itsScf->Iterate(params);
    // Snapshot the (cumulative, process-wide) integrals cache into the run's `cache` section, inside the bracket.
    qchem::BasisSet::theCache<double>().EmitReport();
    RebuildSampling();
    return ok;
}

void AtomCalculation::RebuildSampling()
{
    const auto* wf = itsScf->GetWaveFunction();
    itsDensity = wf->GetChargeDensity();               // it BUILDS the density and hands ownership over

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
    // V1.25: this LEAKED before -- GetChargeDensity BUILDS a whole composite density, and the raw pointer
    // was dropped after reading one number off it.  The owning return makes the temporary free itself.
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

BasisSet::irrepv_t AtomCalculation::GetIrreps(const Spin& ms) const
{
    return itsBasis->GetIrreps(ms);
}

const AtomCalculation::orbital_t* AtomCalculation::GetOrbital(size_t index, const Irrep& qns) const
{
    // Iterate<orbital_t>() yields the concrete orbitals in level order (matches the old QchemTester
    // accessor) -- the plain Iterate() order is not the per-irrep level order the Dirac diagnostics want.
    const auto* orbs = itsScf->GetWaveFunction()->GetOrbitals(qns);
    for (const auto* o : orbs->template Iterate<orbital_t>())
        if (index-- == 0) return o;
    return nullptr;
}

} // namespace qchem
