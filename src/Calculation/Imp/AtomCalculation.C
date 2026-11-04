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
import qchem.BasisSet.IntegralsCacheReport;          // EmitIntegralsCacheReport (the run's `cache` section)
import qchem.Reporting;                              // report:: run sink + the scf-section writer
import qchem.SCFRecipeReport;                        // ReportSCFRecipe: the `scf` header every facade files (V-FACADE-PARITY)

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

// The atom's own keys of the run report's `scf` header (the knobs are the shared ReportSCFRecipe's): what model, which spin
// representation, how many electrons, and a pseudo-atom's valence.
static qchem::report::json AtomRecipeKeys(const AtomCalcOptions& o, int ne)
{
    qchem::report::json own;
    own["model"]      = o.pseudopotential ? std::string("PP") : ModelName(o.model);
    own["pol"]        = (o.spin == SpinGroup::Polarized) ? "P" : "U";
    own["nElectrons"] = ne;
    if (o.pseudopotential) own["valence"] = o.valence;
    return own;
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
        return BasisSet::Radial::Factory(nlohmann::json{{"type", opts.type}, {"shells", shells}}, ec);
    }
    // Explicit single exponent list ("bring your own exponents") -- applied to every occupied l (with ltrim).
    if (!opts.exponents.empty())
        return BasisSet::Radial::Factory(nlohmann::json{{"type", opts.type},
                                                      {"exponents", opts.exponents}, {"ltrim", opts.ltrim}}, ec);
    if (opts.N > 0)
        return BasisSet::Radial::Factory(nlohmann::json{{"type", opts.type}, {"N", opts.N},
                                                      {"emin", opts.emin}, {"emax", opts.emax}}, ec);
    return BasisSet::Radial::Factory(opts.basisPool, opts.type, Z, ec);
}

AtomCalculation::AtomCalculation(int Z, int charge, const AtomCalcOptions& opts,
                                 const SCFParams& params, const AcceleratorOptions& acc)
    : AtomCalculation(BuildOnly{}, Z, charge, opts, acc, nullptr)
{
    Converge(params);
}

std::unique_ptr<AtomCalculation> AtomCalculation::Prepare(int Z, int charge, const AtomCalcOptions& opts,
                                                          const AcceleratorOptions& acc, Observer onIteration)
{
    std::unique_ptr<AtomCalculation> c(new AtomCalculation(BuildOnly{}, Z, charge, opts, acc, std::move(onIteration)));
    c->BuildGraph();                                  // into the CALLER'S run (user ruling D, 2026-11-02: drivers own the runs)
    c->itsPrepared=true;
    return c;
}

AtomCalculation::AtomCalculation(BuildOnly, int Z, int charge, const AtomCalcOptions& opts, const AcceleratorOptions& acc,
                                 Observer onIteration)
    : itsZ(Z)
    , itsNe(Z - charge)
    , itsOpts(opts)
    , itsAcc(acc)
    , itsStructure(std::make_shared<Atom>(Z, double(charge), rvec3_t(0,0,0)))
    , itsObserver(std::move(onIteration))
{
    // The EC follows the model: a pseudo-atom uses the valence-only PseudoAtom_EC (keyed by the full Z and
    // the net charge Zion-electrons); relativistic models need the Dirac (j-coupled) configuration;
    // everything else is a plain Atom_EC.
    itsEC = opts.pseudopotential ? static_cast<ElectronConfiguration*>(new PseudoAtom_EC(itsZ, PPZion(opts,itsNe)-itsNe, PPZion(opts,itsNe)))
          : IsDirac(opts.model)  ? static_cast<ElectronConfiguration*>(new AtomDirac_EC(itsNe))
          :                        static_cast<ElectronConfiguration*>(new Atom_EC(itsNe));
    // The basis is NOT built here: BuildGraph builds it INSIDE an open run report's "basis" scope, so the shells'
    // exponents are announced by the shells themselves as they are constructed (V1.14).
}

AtomCalculation::~AtomCalculation()
{
    delete itsScf;     // FIRST: the iterator references itsHam/itsAccel, which are freed after this body (R2.22)
    delete itsBasis;
    delete itsEC;
}

// The accelerator TYPE the caller asked for (opts.accelerator's "type"; DIIS by default) -- the scf header names it.
std::string AtomCalculation::AcceleratorTag() const
{
    return (itsOpts.accelerator.is_object() && !itsOpts.accelerator.empty())
         ? itsOpts.accelerator.value("type", std::string("DIIS")) : std::string("DIIS");
}

// THE GRAPH, built inside the caller's open report run (if any: with none, every report write is a no-op).
void AtomCalculation::BuildGraph()
{
    namespace H = qchem::Hamiltonian;
    // Atoms use the proven Z-scaled DIIS engage gate unless the caller pins one.  The json escape hatch (opts.accelerator)
    // then merges over these defaults and selects the accelerator TYPE -- the accelerator-tuning surface the scfrun driver
    // needs (DIIS/GDM/Ladder/directmin).  Z^2*0.08/32, was 0.1/32: the [F,D] it is compared against changed UNIT on
    // 2026-10-29 (the worst irrep, was sqrt(sum of squares) over the l-blocks -- user ruling), and on Kr the two differ by
    // ~0.8 at the engage step; unrescaled, DIIS engaged one step early and the Kr Libxc run went chaotic.
    const double emax = itsAcc.eMax > 0.0 ? itsAcc.eMax : itsZ*itsZ*0.08/32;
    nlohmann::json jsacc = ToAcceleratorJson(itsAcc, emax);   // every knob of the chain; the escape hatch below merges over it
    const std::string ts = AcceleratorTag();
    if (itsOpts.accelerator.is_object() && !itsOpts.accelerator.empty())
        for (auto& [k,v] : itsOpts.accelerator.items()) jsacc[k]=v;
    using AType = qchem::SCFAccelerators::Type;
    const bool directmin = (ts=="directmin" || ts=="DirectMin");
    if (directmin) jsacc["EMax"]=1e10;   // GDM always steps; the direct-min loop seeds via diagonalize
    const AType atype = directmin ? AType::GDM : ts=="Ladder" ? AType::Ladder : ts=="GDM" ? AType::GDM : AType::DIIS;
    // The facade sets the context, the accelerator files what it was built with (THE RUN REPORT): one `accelerators` row per
    // accelerator built in this run (a solid rebuilds one per stage), its Ladder rungs one level further in.
    qchem::SCFAccelerators::SCFAccelerator* accel=nullptr;
    { qchem::report::Row built("accelerators"); accel = qchem::SCFAccelerators::Factory(atype, jsacc); }
    // Default seed for atoms is the core guess (atoms never use the molecular SAD seed).
    using qchem::ChargeDensity::SeedStrategy;
    const auto seed = (itsOpts.seed != SeedStrategy::Default) ? itsOpts.seed : SeedStrategy::CoreGuess;
    // RELATIVISTIC (Dirac) bases are even-tempered / kinetic-balanced: their overlap is near-singular BY
    // CONSTRUCTION, and canonical-ortho mode-dropping breaks kinetic balance -> variational collapse.  So a
    // Dirac SCF ALWAYS uses PLAIN Cholesky (never truncate); a non-relativistic atom uses opts.ortho (Auto by
    // default, but Cholesky when a caller deliberately wants an over-complete basis kept whole).
    const qchem::Ortho ortho = IsDirac(itsOpts.model) ? qchem::Cholesky : itsOpts.ortho;

    namespace rpt = qchem::report;

    // The `basis` section is ASSEMBLED across layers, and since V1.14 it spans the basis's whole LIFE in this
    // run: the FIRST Converge BUILDS the basis (each exponent shell announces its own basis.exponents row as it is
    // constructed -- nobody asks it to), the orchestrator stamps the scalars, and the WF build (inside `new SCFIter`)
    // fills basis.perIrrep/removed via the cursor (the shared LASolver conditioning).  The HAMILTONIAN is built OUTSIDE
    // it, in its own context "grids" (user ruling 2026-11-02: the caller sets the context, the callee reports one level
    // in -- its meshes file there, the same paths as a molecule's and a solid's).  So the basis context is entered twice:
    // SILENTLY for the build (a Scope renders only the scalars it added; the build adds rows), then as the rendering
    // Section, which renders the whole basis once, complete, on close.  A later Converge reuses the basis and truthfully
    // carries no exponents rows.
    {
        rpt::Scope basis("basis");
        if (!itsBasis) itsBasis = BuildBasis(itsOpts, itsZ, *itsEC);
    }
    // Three DFT/HF routes: a pseudopotential (PP front door, valence electrons = itsNe), an explicit XC functional
    // override (the public selector, e.g. libxc), or the model's built-in Hamiltonian/functional.
    decltype(itsHam.get()) ham=nullptr;
    {
        rpt::Scope grids("grids");
        ham = itsOpts.pseudopotential
            ? H::Factory(itsOpts.spin, itsStructure, thePeriodicTable().GetSymbol(itsZ), PPZion(itsOpts,itsNe), itsOpts.mesh, itsBasis,
                         itsOpts.vxcFit)
            : itsOpts.xc.has_value()
                ? H::Factory(itsOpts.spin, itsStructure, *itsOpts.xc, itsOpts.mesh, itsBasis, itsOpts.vxcFit)
                : H::Factory(itsOpts.model, itsOpts.spin, itsStructure, itsOpts.mesh, itsBasis, itsOpts.xalpha, itsOpts.vxcFit);
    }
    {
        rpt::Section basis("basis");
        rpt::Set("type",       AtomTypeName(itsOpts.type));
        rpt::Set("nFunctions", (long)itsBasis->GetNumFunctions());
        // R2.22: the iterator no longer deletes what it is handed, so this facade adopts the pair.  Order is
        // deliberate -- the PREVIOUS iterator dies first, then the previous Hamiltonian/accelerator it pointed
        // at are freed by these resets, which is exactly when `delete itsScf` used to free them.
        delete itsScf; itsScf = nullptr;
        itsHam.reset(ham);
        itsAccel.reset(accel);
        itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), seed, itsStructure.get(), ortho);
    }
    // (directmin => AType::GDM above, whose WantsLineSearch() drives the direct-min loop -- no SetDirectMin needed.)
}

bool AtomCalculation::Converge(const SCFParams& params)
{
    // The atomic SCF writes into the CALLER'S report run, as every facade does (user ruling D, 2026-11-02): scfrun, valgen, the
    // valence-basis generator and the tests that read the report open it.  With no run open every write is a no-op.
    namespace rpt = qchem::report;
    { rpt::Section scf("scf"); ReportSCFRecipe(params, AcceleratorTag(), AtomRecipeKeys(itsOpts, itsNe)); ReportAccuracy(itsOpts.accuracy, kNoAccuracyConsumerYet); }   // the "what recipe am I running" header
    if (itsPrepared) itsPrepared=false;                // Prepare built the graph (in its caller's run): reuse it, once
    else             BuildGraph();
    if (itsObserver) itsScf->SetObserver(itsObserver);

    bool ok = itsScf->Iterate(params);
    // Snapshot the (cumulative, process-wide) integrals cache into the run's `cache` section, inside the bracket.
    qchem::BasisSet::EmitIntegralsCacheReport();
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
