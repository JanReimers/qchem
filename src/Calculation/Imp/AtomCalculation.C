// File: Calculation/Imp/AtomCalculation.C  Implementation of the qchem::AtomCalculation facade.
module;
#include <memory>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include <map>
#include <optional>
#include <iomanip>
module qchem.AtomCalculation;

import qchem.SetupDecisions;   // DIISEngageGate: the engage gate every setup shares (V-FACADE-SETUP)
import qchem.Hamiltonian.Factory;                   // Factory(...), IsDFT, XCFunctional (XC selector)
import qchem.Policy;                                  // policy::PPFittedFunctional (the PP / functional match)
import qchem.Pseudopotential.GTH_Potentials;          // GetGTH(...).fitted: the functional the pseudo-atom's PP was fitted with
import qchem.XCRhoPolicy;                            // ReportXCRho: the XC ρ-axis policy row
import qchem.ElectronConfiguration.AtomNR;          // Atom_EC, PseudoAtom_EC
import qchem.ElectronConfiguration.AtomDirac;       // AtomDirac_EC
import qchem.PeriodicTable;                          // thePeriodicTable().GetSymbol(Z) (the PP element)
import qchem.SCFAccelerator.Factory;                // SCFAccelerators::Type, Factory
import qchem.WaveFunction;                           // WaveFunction (GetChargeDensity/GetOrbitals/GetQNs)
import qchem.Orbitals;                               // Orbital, Orbitals
import qchem.Types;                                  // Vector3D / rvec3_t
import qchem.BasisSet.Conditioning;                  // ReportConditioning
import qchem.Reporting;                              // report:: run sink + the scf-section writer
import qchem.SCFRecipeReport;                        // ReportSCFRecipe: the `scf` header every facade files (V-FACADE-PARITY)
import qchem.BasisSet.Orbital_1E_IBS;                // Real_OIBS -- the per-irrep blocks a saved state walks
import qchem.CompositeCD;                            // tComposite_CD (a restored density's shape)
import qchem.ChargeDensity.Factory;                  // IrrepCD_Factory
import qchem.BasisSet.Radial.Factory;                // CoulombFitLadder: the Dunlap fit's rungs (the pilot's third axis)
import qchem.Mesh.BeckeClimb;                        // BeckeTailEstimate: a rung's error from its difference to the rung above
import qchem.Accuracy;                               // AccuracyTarget, ClientConsumes
import qchem.ERI4MemoryPolicy;                       // policy::ERI4Memory: the ERI4 RAM cap at setup (S6)
import qchem.Math;                                   // std::fabs, std::max (the pilot's distances)

namespace qchem
{

using SCFIter = qchem::SCFIterator::MolecularSCFIterator;   // item 2: the atom/molecule display subclass

static bool IsDirac(Model m) {return m == Model::DE1 || m == Model::DHF;}

// Human-readable model tag for the run report's `scf` section (a PP run reports "PP" -- it bypasses the Model enum).
static std::string ModelName(Model m) {return qchem::Hamiltonian::Name(m);}

//! A pseudo-atom's functional: an \c xc override when stated, PADE when the model says so (the GTH-PADE PPs' own, 2026-11-12), else the
//! historical LDA (Dirac + VWN5) -- a PP atom has always ignored any other model, and still does (its tests leave model at HF).
static qchem::Hamiltonian::XCFunctional PseudoAtomXC(const AtomCalcOptions& o)
{
    if (o.xc) return *o.xc;
    return o.model==Model::PADE ? qchem::Hamiltonian::ModelToXC(Model::PADE) : qchem::Hamiltonian::XCFunctional{};
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

// THE CONFIGURATION POLICY (OpenWork row "All-electron ION seeds", 2026-11-15): STATED (opts.configuration, NIST's notation) or derived --
// the isoelectronic neutral all-electron, the pseudo-ion rule under a pseudopotential, the j-coupled one for Dirac.  Either way the run
// report files the configuration the SCF occupies (policy AtomConfiguration: Choose when derived; a stated one that cannot run is REFUSED).
ElectronConfiguration* AtomCalculation::MakeConfiguration() const
{
    const std::string el=thePeriodicTable().GetSymbol(itsZ);
    const int zion=PPZion(itsOpts,itsNe);
    const int core=itsOpts.pseudopotential ? itsZ-zion : 0;   // the electrons the SCF does not carry
    if (itsOpts.configuration.empty())
    {
        if (IsDirac(itsOpts.model)) return new AtomDirac_EC(itsNe);
        Atom_EC* ec=itsOpts.pseudopotential ? new PseudoAtom_EC(itsZ, zion-itsNe, zion) : new Atom_EC(itsNe);
        policy::Report("AtomConfiguration", el+" "+ConfigurationText(ec->Shells()), policy::Verdict::Choose,
                       itsOpts.pseudopotential ? "unstated: the pseudo-ion rule (a cation empties s before d)"
                                               : "unstated: the isoelectronic neutral's configuration (right for an anion and a main-group cation; "
                                                 "a d- or f-block cation states its own)");
        return ec;
    }
    const std::string stated="configuration \""+itsOpts.configuration+"\"";
    if (IsDirac(itsOpts.model))
        policy::Refuse("AtomConfiguration", stated+" on a Dirac atom", "a stated configuration on a j-coupled (Dirac) atom is not realized yet "
                       "(doc/OpenWork.md); leave it unstated for the derived one");
    const auto parsed=ParseConfiguration(itsOpts.configuration);
    if (!parsed) policy::Refuse("AtomConfiguration", stated, parsed.Error());
    const StatedConfiguration& c=parsed.Value();
    // A pseudo-atom's [X] is the core its pseudopotential REPLACES (the SCF carries only the named shells); an all-electron atom carries its
    // [X] too, so there the electron count below is the whole check.
    if (itsOpts.pseudopotential && c.CoreElectrons()!=core)
        policy::Refuse("AtomConfiguration", stated+" for "+el+" q"+std::to_string(zion),
                       "its ["+c.coreSymbol+"] core holds "+std::to_string(c.CoreElectrons())+" electrons, but the pseudopotential replaces "
                       +std::to_string(core)+" (Z - q): write the core it replaces and name the rest");
    const ShellOccupations carried=itsOpts.pseudopotential ? c.named : c.All();
    int n=0; for (const auto& o : carried) n+=o.count;
    if (n!=itsNe)
        policy::Refuse("AtomConfiguration", stated+" for "+el+" at charge "+std::to_string(itsZ-itsNe-core),
                       "it holds "+std::to_string(n)+" electrons where the run carries "+std::to_string(itsNe));
    if (const std::string why=WhyNotRepresentable(carried, /*fromLowestN*/!itsOpts.pseudopotential); !why.empty())
        policy::Refuse("AtomConfiguration", stated, why);
    qchem::report::Set("configuration", ConfigurationText(carried));   // a stated one is the caller's: no policy row, but the record carries it
    return new Atom_EC(carried);
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
    itsEC = MakeConfiguration();
    itsDFT = opts.pseudopotential || opts.xc.has_value() || qchem::Hamiltonian::IsDFT(opts.model);
    // THE RADIAL-ONLY ALIAS (V-XC-ONE-TERM step 5, controller 2026-11-14: deck backward compatibility).  Before AngularKind::None the atom
    // spelled its radial-only set Lebedev at degree 0, the 1-direction rule: the SAME points (one direction at 4 pi).  On the atom that
    // spelling still reads as it always did, now under its own name, so the admission rule and the reports see what the mesh IS.
    if (itsOpts.mesh.angular==qcMesh::AngularKind::Lebedev && itsOpts.mesh.angularDegree==0) itsOpts.mesh.angular=qcMesh::AngularKind::None;
    if (itsOpts.mesh.angular==qcMesh::AngularKind::None) itsOpts.mesh.angularDegree=0;   // not read by the rule: one value, one mesh ID
    ResolveGrids();
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
    const double emax = DIISEngageGate(itsAcc.eMax, itsZ, /*the atom keeps the Z-scaled gate, DFT too*/false);
    nlohmann::json jsacc = ToAcceleratorJson(itsAcc, emax);   // every knob of the chain; the escape hatch below merges over it
    const std::string ts = AcceleratorTag();
    if (itsOpts.accelerator.is_object() && !itsOpts.accelerator.empty())
        for (auto& [k,v] : itsOpts.accelerator.items()) jsacc[k]=v;
    if (itsStartOnGDM) { jsacc["startOnGDM"]=true; itsStartOnGDM=false; }   // the restored stage only (PrepareRestart)
    using AType = qchem::SCFAccelerators::Type;
    const bool directmin = (ts=="directmin" || ts=="DirectMin");
    if (directmin) jsacc["EMax"]=1e10;   // GDM always steps; the direct-min loop seeds via diagonalize
    const AType atype = directmin ? AType::GDM : ts=="Ladder" ? AType::Ladder : ts=="GDM" ? AType::GDM : AType::DIIS;
    if (itsOpts.pseudopotential || itsOpts.xc.has_value() || H::IsDFT(itsOpts.model))   // GDM needs H = dE/dD (pin 32)
        RefuseGradientConsumers(itsOpts.vxcFit, itsOpts.xcRho, atype==AType::GDM ? "GDM" : atype==AType::Ladder ? "Ladder" : "DIIS");
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
    // constructed -- nobody asks it to), the orchestrator stamps the scalars and asks for the conditioning
    // report (basis.perIrrep, removed, conditioning).  The HAMILTONIAN is built OUTSIDE
    // it, in its own context "hamiltonian" (user rulings 2026-11-02, 2026-11-04: the caller sets the context, the callee
    // reports one level in -- its meshes and terms file there, the same paths as a molecule's and a solid's).  So the basis context is entered twice:
    // SILENTLY for the build (a Scope renders only the scalars it added; the build adds rows), then as the rendering
    // Section, which renders the whole basis once, complete, on close.  A later Converge reuses the basis and truthfully
    // carries no exponents rows.
    {
        rpt::Scope basis("basis");
        if (!itsBasis) itsBasis = BuildBasis(itsOpts, itsZ, *itsEC);
    }
    // THE ERI4 RAM CAP (S6, user ruling 6): an all-electron HF / DHF atom builds the 4-centre tensors; judged once here from the block sizes.
    // (A PP atom ignores the model and runs LDA, and an XC override is DFT: neither builds ERI4.)
    if (!itsOpts.pseudopotential && !itsOpts.xc.has_value() && (itsOpts.model==Model::HF || itsOpts.model==Model::DHF))
    {
        std::vector<size_t> blocks;
        for (auto b : itsBasis->Iterate<BasisSet::Real_IBS>()) blocks.push_back(b->GetNumFunctions());
        policy::ERI4Memory(blocks, "atom Z="+std::to_string(itsZ));
    }
    // Three DFT/HF routes: a pseudopotential (PP front door, valence electrons = itsNe), an explicit XC functional
    // override (the public selector, e.g. libxc), or the model's built-in Hamiltonian/functional.
    decltype(itsHam.get()) ham=nullptr;
    {
        rpt::Scope hamiltonian("hamiltonian");
        if (itsOpts.pseudopotential || itsOpts.xc.has_value() || H::IsDFT(itsOpts.model))
            ReportXCRho(itsOpts.vxcFit, itsOpts.xcRho);   // the ρ axis (user ruling 2026-11-09): warns on the Gaussian x δ cell
        if (itsOpts.pseudopotential)   // THE PP / FUNCTIONAL MATCH (user ruling 2026-11-12): warn and run, never switch
        {
            const std::string el=thePeriodicTable().GetSymbol(itsZ);
            const auto xcRun=PseudoAtomXC(itsOpts);
            policy::PPFittedFunctional(H::PPFunctionalChecks({{el, Pseudopotential::GetGTH(el, "LDA", PPZion(itsOpts,itsNe)).fitted}}, xcRun));
        }
        ham = MakeHamiltonianOn(itsGrids, *itsCoulombFit);
    }
    {
        rpt::Section basis("basis");
        rpt::Set("type",       AtomTypeName(itsOpts.type));
        rpt::Set("nFunctions", (long)itsBasis->GetNumFunctions());
        qchem::BasisSet::ReportConditioning(*itsBasis);
        // R2.22: the iterator no longer deletes what it is handed, so this facade adopts the pair.  Order is
        // deliberate -- the PREVIOUS iterator dies first, then the previous Hamiltonian/accelerator it pointed
        // at are freed by these resets, which is exactly when `delete itsScf` used to free them.
        delete itsScf; itsScf = nullptr;
        itsHam.reset(ham);
        itsAccel.reset(accel);
        if (itsCarried)   // the atom pilot's re-convergence on raised grids: the converged D drives the first Fock
            itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), itsCarried.release(), itsStructure.get(), ortho,
                                 0.0, qchem::SCFIterator::SeedRole::WorkingDensity);
        else if (itsRestored.empty())
            itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), seed, itsStructure.get(), itsOpts.mesh,
                                 // a stated SAD seed reads the run's own atom: the GTH pseudo-atom at its Zion (the valence library), else all-electron
                                 itsOpts.pseudopotential ? qchem::ChargeDensity::SeedKey("LDA", {{itsZ, PPZion(itsOpts,itsNe)}}) : qchem::ChargeDensity::SeedKey{},
                                 ortho);   // an HF SAD seed's sibling: the run's mesh
        else
        {
            // A RESTORED state (PrepareRestart), consumed here once: an exact resume drives the first Fock from the saved D.
            itsScf = new SCFIter(itsBasis, itsEC, itsHam.get(), itsAccel.get(), RestoredDensity().release(), itsStructure.get(), ortho,
                                 0.0, itsRestoredExact ? qchem::SCFIterator::SeedRole::ExactResume : qchem::SCFIterator::SeedRole::WorkingDensity);
            itsRestored.clear();
            itsRestoredExact=false;
        }
    }
    // (directmin => AType::GDM above, whose WantsLineSearch() drives the direct-min loop -- no SetDirectMin needed.)
}

bool AtomCalculation::Converge(const SCFParams& params)
{
    // The atomic SCF writes into the CALLER'S report run, as every facade does (user ruling D, 2026-11-02): scfrun, valgen, the
    // valence-basis generator and the tests that read the report open it.  With no run open every write is a no-op.
    namespace rpt = qchem::report;
    { rpt::Section scf("scf"); ReportSCFRecipe(params, AcceleratorTag(), AtomRecipeKeys(itsOpts, itsNe)); ReportAccuracy(itsOpts.accuracy, itsDFT && itsOpts.accuracy.eps>0.0
            ? std::string("the XC radial mesh and the Hartree route's mesh or Coulomb fit (the atom pilot, at the converged density)")
            : std::string(kNoAccuracyConsumerYet)); }   // the "what recipe am I running" header
    if (itsPrepared) itsPrepared=false;                // Prepare built the graph (in its caller's run): reuse it, once
    else             BuildGraph();
    if (itsObserver) itsScf->SetObserver(itsObserver);

    itsLastkT=params.SmearingkT;
    bool ok = itsScf->Iterate(params);
    if (ok && itsDFT && !itsPiloted && itsOpts.accuracy.eps>0.0) { RebuildSampling(); PilotGrids(params); ok=itsScf->Converged(); }
    // (No end-of-run cache snapshot since S6: each owner-held table reported its own bytes under `cache` where it was filled.)
    RebuildSampling();
    ReportResult("atom "+thePeriodicTable().GetSymbol(itsZ), ok, itsScf->GetIterationCount(), itsDensity->GetTotalCharge(), itsScf->GetEnergy());
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

//====================================================================================================
//  THE SAVED STATE (V-FACADE-STATE): the atom's hooks and its restart
//====================================================================================================
Calculation::StateIdentity AtomCalculation::MakeStateIdentity() const
{
    return {.Nelec=itsNe, .spin=itsOpts.spin, .multiplicity=0, .functional=FunctionalName(), .kT=itsLastkT};
}

std::string AtomCalculation::FunctionalName() const
{
    if (itsOpts.pseudopotential)   // the PP front door: LDA (Dirac + VWN5), or PADE when stated
        return itsOpts.xc ? "LDA (xc override)" : itsOpts.model==Model::PADE ? "PADE" : "LDA";
    return ModelName(itsOpts.model)+(itsOpts.xc ? " (xc override)" : "");
}

std::vector<StateBlock> AtomCalculation::StateBlocks(bool withOrbitals) const
{
    if (!itsBasis) throw std::logic_error("AtomCalculation: no basis yet (a saved state needs the graph)");
    const auto* wf = withOrbitals && itsScf ? itsScf->GetWaveFunction() : nullptr;
    if (withOrbitals && !wf) throw std::logic_error("AtomCalculation: no wave function yet (a saved state needs a Converge)");
    std::vector<StateBlock> out;   // spin irreps outer, basis blocks inner (the order the composite wave function builds)
    for (Spin s : SpinIrreps(itsOpts.spin))
        for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            StateBlock b;
            b.real=ibs;
            b.irrep=ibs->GetIrrep(s);
            if (wf) b.orbitals=wf->GetOrbitals(b.irrep);
            out.push_back(std::move(b));
        }
    return out;
}

StateSummary AtomCalculation::MakeStateSummary() const
{
    if (!itsScf || !itsScf->GetWaveFunction()) throw std::logic_error("AtomCalculation: no wave function yet");
    return {.label="atom "+thePeriodicTable().GetSymbol(itsZ), .converged=itsScf->Converged(), .iterations=itsScf->GetIterationCount(),
            .energy=itsScf->GetEnergy().GetTotalEnergy(), .charge=itsDensity ? itsDensity->GetTotalCharge() : 0.0,
            .recipe="model "+FunctionalName()+", basis "+AtomTypeName(itsOpts.type)+", accelerator "+AcceleratorTag()};
}

void AtomCalculation::RecordState(StateRecord& r, double kT) const
{
    Calculation::RecordState(r, kT);   // the structure, electron and block keys, and the blocks
    const AtomCalcOptions& o=itsOpts;
    StateFingerprint& fp=r.fp;
    // ---- REFUSE: the saved D is a matrix over THIS radial basis's functions (the blocks' sizes are the shared keys) ----
    fp.refuseStr["basis.family"]="radial "+AtomTypeName(o.type);
    fp.refuseNum["basis.pool"]={double(int(o.basisPool)), double(o.N), o.emin, o.emax, double(o.ltrim)};
    fp.refuseNum["basis.exponents"]=o.exponents;
    {
        std::vector<double> byL;   // l, count, exponents... per entry
        for (const auto& [l,es] : o.exponentsByL) { byL.push_back(l); byL.push_back(es.size()); byL.insert(byL.end(), es.begin(), es.end()); }
        fp.refuseNum["basis.exponentsByL"]=byL;
    }
    fp.refuseNum["pp.valence"]={o.pseudopotential ? double(PPZion(o, itsNe)) : 0.0};
    if (!o.configuration.empty()) fp.refuseStr["configuration"]=o.configuration;   // a state of another configuration seeds nothing (stated only: the derived default keeps its old fingerprint)
    // ---- WARM ----
    if (o.model==Model::Xalpha) fp.warmNum["xc.xalpha"]={o.xalpha};
    fp.warmNum["grids.xcMesh"]={double(o.mesh.nRadial), double(o.mesh.angularDegree), double(int(o.mesh.radial)), double(int(o.mesh.angular))};
    fp.warmNum["xc.fit"]={double(int(o.vxcFit))};
    fp.warmNum["coulomb.route"]={double(int(o.coulomb.fit)), double(int(o.coulomb.solver)), double(o.coulomb.lMax), double(o.coulomb.nRadial)};
    fp.warmNum["ortho.kind"]={double(int(o.ortho))};
}

const std::optional<RestartSummary>& AtomCalculation::RestartedFrom() const {return itsRestartedFrom;}

std::unique_ptr<qchem::ChargeDensity::rDM_CD> AtomCalculation::RestoredDensity() const
{
    using namespace qchem::ChargeDensity;
    auto cd=std::make_unique<tComposite_CD<double>>(itsBasis->GetImposition() ? itsBasis->GetImposition()->DensityGroup() : nullptr);
    size_t b=0;
    for (Spin s : SpinIrreps(itsOpts.spin))
        for (auto ibs : itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            if (b>=itsRestored.size() || itsRestored[b].rows()!=ibs->GetNumFunctions())
                throw std::logic_error("AtomCalculation: the restored density's blocks do not match this basis");
            const Irrep irr=ibs->GetIrrep(s);
            cd->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(itsRestored[b++], ibs, irr)), irr);
        }
    if (b!=itsRestored.size()) throw std::logic_error("AtomCalculation: the restored density has more blocks than this basis");
    return cd;
}

Outcome<std::unique_ptr<StateView>, RestartRefusal>
AtomCalculation::ViewState(const std::string& path, int Z, int charge, const AtomCalcOptions& opts, const AcceleratorOptions& acc)
{
    using O=Outcome<std::unique_ptr<StateView>, RestartRefusal>;
    auto saved=ReadState(path);
    if (!saved) return O::Fail(saved.Error());
    std::unique_ptr<AtomCalculation> c(new AtomCalculation(BuildOnly{}, Z, charge, opts, acc, nullptr));
    {
        qchem::report::Scope basis("basis");
        c->itsBasis=BuildBasis(c->itsOpts, c->itsZ, *c->itsEC);
    }
    auto judged=c->JudgeView(*saved);
    if (!judged) return O::Fail(judged.Error());
    std::vector<StateBlock> blocks=c->StateBlocks(false);
    auto D=CheckedBlockDensities(*saved, blocks);
    if (!D) return O::Fail(D.Error());
    for (const mat_t<dcmplx>& M : *D)
    {
        const size_t n=M.rows();
        hmat_t<double> h(n);
        for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) h(a,b)=M(a,b).real();
        c->itsRestored.push_back(std::move(h));
    }
    std::unique_ptr<const qchem::ChargeDensity::rDM_CD> rho=c->RestoredDensity();
    const auto* up  =qchem::ChargeDensity::ChannelOf<double>(rho.get(), Spin::Up);
    const auto* down=qchem::ChargeDensity::ChannelOf<double>(rho.get(), Spin::Down);
    if (!up || !down) up=down=nullptr;   // an unpolarized state resolves no channel: m is the zero field
    return O::Ok(std::make_unique<StateView>(std::move(c), saved.TakeValue(), judged.TakeValue(), std::move(blocks), std::move(rho), up, down));
}

Outcome<std::unique_ptr<AtomCalculation>, RestartRefusal>
AtomCalculation::Restart(const std::string& path, int Z, int charge, const AtomCalcOptions& opts, const SCFParams& params,
                         const AcceleratorOptions& acc, Observer onIteration)
{
    auto R=PrepareRestart(path, Z, charge, opts, params, acc, std::move(onIteration));
    if (R) (void)R.Value()->Converge(params);
    return R;
}

Outcome<std::unique_ptr<AtomCalculation>, RestartRefusal>
AtomCalculation::PrepareRestart(const std::string& path, int Z, int charge, const AtomCalcOptions& opts, const SCFParams& params,
                                const AcceleratorOptions& acc, Observer onIteration)
{
    using O=Outcome<std::unique_ptr<AtomCalculation>, RestartRefusal>;
    const std::string tag="[atom restart] ";
    auto saved=ReadState(path);
    if (!saved)
    {
        std::cout<<tag<<"REFUSED ("<<path<<"): "<<saved.Error().details<<std::endl;
        return O::Fail(saved.Error());
    }
    std::unique_ptr<AtomCalculation> c(new AtomCalculation(BuildOnly{}, Z, charge, opts, acc, std::move(onIteration)));
    {
        qchem::report::Scope basis("basis");
        c->itsBasis=BuildBasis(c->itsOpts, c->itsZ, *c->itsEC);   // BuildGraph reuses it
    }
    SavedState state=saved.TakeValue();
    auto judged=c->JudgeRestart(state, params.SmearingkT, tag);
    if (!judged) return O::Fail(judged.Error());
    auto D=CheckedBlockDensities(state, c->StateBlocks(false));
    if (!D)
    {
        std::cout<<tag<<"REFUSED ("<<path<<"): "<<D.Error().details<<std::endl;
        return O::Fail(D.Error());
    }
    for (const mat_t<dcmplx>& M : *D)
    {
        const size_t n=M.rows();
        hmat_t<double> h(n);
        for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) h(a,b)=M(a,b).real();
        c->itsRestored.push_back(std::move(h));
    }
    c->itsRestoredExact=judged->exact;
    c->itsStartOnGDM=judged->exact && params.SmearingkT==0.0 && c->AcceleratorTag()=="Ladder";
    c->itsRestartedFrom=*judged;
    c->BuildGraph();
    c->itsPrepared=true;
    return O::Ok(std::move(c));
}

} // namespace qchem

namespace qchem
{
// ====================================================================================================================================
// THE ATOM'S GRIDS (pin 31; lane B 2026-11-14, the atom-below-NIST row).  Two approximations were fixed and eps-blind, and both pulled
// the energy DOWN: the 50-point XC radial mesh (Xe 0.6-1.7 mHa, U 1-6) and the Dunlap Coulomb fit on the s block x 2 (its robust energy
// is J - (Δρ|Δρ)/2 <= J: Slater U 32 mHa).  Converged, every atom lands above NIST.  So both are SIZED here, from the accuracy target.
// ====================================================================================================================================
namespace
{
//! The ladders the pilot climbs (policy, measured on Xe and U, Slater and Gaussian High, 2026-11-14): the XC count is flat from 100 on
//! Xe and from 150 on U Gaussian; the Hartree's radial-Poisson mesh moves U 40 μHa from 200 to 400 and nothing from 400 to 800.
constexpr int kAtomXCRadials[]     ={50, 75, 100, 150, 200, 300, 400};
constexpr int kAtomHartreeRadials[]={200, 300, 400, 600, 800, 1200};
constexpr int kAtomPilotRounds=3;   //!< the first round climbs freely; each later one confirms at the re-converged D, raising at most once
template <size_t N> std::optional<int> RungAbove(const int (&rungs)[N], int n)
{
    for (int r : rungs) if (r>n) return r;
    return std::nullopt;
}
}

void AtomCalculation::ResolveGrids()
{
    namespace H = qchem::Hamiltonian;
    if (!itsDFT) return;   // HF, Dirac, 1-electron: no mesh, no Hartree route
    const bool sized=itsOpts.accuracy.eps>0.0;
    itsGrids.xcStated=itsOpts.mesh.nRadial>0;
    itsGrids.xcRadial=itsGrids.xcStated ? itsOpts.mesh.nRadial : kAtomXCRadials[0];
    itsGrids.coulomb=itsOpts.coulomb;
    using CF=H::CoulombFit;
    if (itsGrids.coulomb.fit==CF::Delta)
        policy::Refuse("AtomHartreeRoute", "the δ Hartree route on an atom",
                       "its radial-Poisson projection needs angular channels, and the atom's mesh is the radial-only set (one direction), so every channel "
                       "above l = 0 aliases: Xe gives half its energy and never converges.  The Multipole route is the δ route's exact l = 0 form "
                       "on this rule -- state Multipole or Auxiliary, or leave the route unstated (Multipole)");
    if (itsGrids.coulomb.fit==CF::Auto)
    {
        itsGrids.coulomb.fit=CF::Multipole;
        policy::Report("AtomHartreeRoute", "the atom's Hartree is the radial-Poisson multipole route (l = 0)", policy::Verdict::Choose,
                       "no Hartree route stated: radial Poisson on the density's l = 0 channel is exact to its mesh and as cheap as the fit "
                       "(U 0.27 s vs 0.36 s); the Dunlap fit (Auxiliary) is one stated choice away");
    }
    if (itsGrids.coulomb.fit==CF::Multipole)
    {
        if (const auto a=H::AdmitCoulomb(*itsStructure, itsGrids.coulomb, itsOpts.mesh); !a)
            policy::Refuse("AdmitCoulomb", "the stated Hartree route on this atom", a.Error());
        itsGrids.hStated=itsGrids.coulomb.nRadial>0;
        if (!itsGrids.hStated) itsGrids.coulomb.nRadial=kAtomHartreeRadials[0];
    }
    itsGrids.fitRung=0;
    const std::string xc=std::to_string(itsGrids.xcRadial);
    if (!sized)
        policy::Report("AtomGridsAsStated", "the atom's grids as stated: XC radial "+xc, policy::Verdict::Choose,
                       "accuracy eps = 0 states no target, so nothing is sized: the XC mesh enters at its rung, and the Dunlap fit (if stated) "
                       "is the s block x 2 (rung 0)");
    else if (!itsGrids.xcStated)
        policy::Report("AtomXCRadialClimbs", "the atom's XC radial mesh climbs from "+xc+" points to eps at the converged density",
                       policy::Verdict::Choose, "no XC radial count stated: the run sizes it from the accuracy target (pin 31), only up");
}

qchem::Hamiltonian::rHamiltonian* AtomCalculation::MakeHamiltonianOn(const Grids& g, const Hamiltonian::CoulombFitSlot<double>& slot) const
{
    namespace H = qchem::Hamiltonian;
    qcMesh::MeshParams mesh=itsOpts.mesh;
    mesh.nRadial=g.xcRadial;
    if (itsDFT && g.coulomb.fit==H::CoulombFit::Auxiliary && g.fitRung>0)
    {
        // The Dunlap fit above rung 0 (rung 0 is the basis's own CreateCDFitBasisSet, which the Hartree term fills itself).  A cross-cast
        // between two abstract faces: the radial basis realizes the ladder; a basis without one never reached a raised rung (PilotGrids).
        const auto* ladder=dynamic_cast<const BasisSet::Radial::CoulombFitLadder*>(itsBasis);
        if (!ladder) throw std::logic_error("AtomCalculation::MakeHamiltonianOn: a raised Coulomb-fit rung on a basis without a CoulombFitLadder");
        using fit_t=std::shared_ptr<const BasisSet::rFIT_CD_ABS>;
        slot.GetOrFill([&]{return fit_t(ladder->MakeCoulombFitRung(g.fitRung));});
    }
    return itsOpts.pseudopotential
        ? H::Factory(itsOpts.spin, itsStructure, thePeriodicTable().GetSymbol(itsZ), PPZion(itsOpts,itsNe), mesh, itsBasis,
                     slot, H::XCRoute(itsOpts.vxcFit, itsOpts.xcRho), g.coulomb, PseudoAtomXC(itsOpts))
        : itsOpts.xc.has_value()
            ? H::Factory(itsOpts.spin, itsStructure, *itsOpts.xc, mesh, itsBasis, slot, H::XCRoute(itsOpts.vxcFit, itsOpts.xcRho), g.coulomb)
            : H::Factory(itsOpts.model, itsOpts.spin, itsStructure, mesh, itsBasis, slot, itsOpts.xalpha,
                         H::XCRoute(itsOpts.vxcFit, itsOpts.xcRho), g.coulomb);
}

// THE ATOM PILOT (pin 31; the molecule's mesh pilot, one-dimensional per axis).  The run has converged on its entry grids.  At THIS D,
// each unstated axis climbs its ladder: the distance between a rung and the rung above is the worst change of any term energy (each axis
// moves only its own term), and the rung's error is BeckeTailEstimate of it.  Up to three axes: the XC radial count, the Hartree route's
// radial count (Multipole) or the Dunlap fit's rung (Auxiliary).  A stated axis is audited once and never raised.  If anything was raised,
// rebuild on the raised grids, re-converge ONCE from the carried D, and confirm (each later round raises an axis at most once).
void AtomCalculation::PilotGrids(const SCFParams& params)
{
    namespace rpt = qchem::report;
    namespace H = qchem::Hamiltonian;
    itsPiloted=true;
    const double target=itsOpts.accuracy.Target(ClientConsumes::Energy, 1);
    const auto* ladder=dynamic_cast<const BasisSet::Radial::CoulombFitLadder*>(itsBasis);
    const bool  dunlap=itsGrids.coulomb.fit==H::CoulombFit::Auxiliary;
    const std::string tag="[atom pilot] ";
    report::FormatGuard keep(std::cout);
    rpt::Section pilot("atomPilot");
    rpt::Set("measure", "E"); rpt::Set("target", target);
    if (dunlap && !ladder)
        policy::Report("AtomCoulombFitNotSized", "the Dunlap fit stays at the basis's own (rung 0)", policy::Verdict::Warn,
                       "this radial basis carries no CoulombFitLadder, so its fit cannot be climbed; its energy error is -(Δρ|Δρ)/2, unmeasured");
    enum class Axis {XC, Hartree, Fit};
    auto name=[](Axis a){return a==Axis::XC ? std::string("xcRadial") : a==Axis::Hartree ? std::string("hartreeRadial") : std::string("coulombFitRung");};
    auto value=[](const Grids& g, Axis a){return a==Axis::XC ? g.xcRadial : a==Axis::Hartree ? g.coulomb.nRadial : g.fitRung;};
    auto above=[&](const Grids& g, Axis a) -> std::optional<Grids>
    {
        Grids u=g;
        if (a==Axis::XC)      { auto n=RungAbove(kAtomXCRadials, g.xcRadial);           if (!n) return std::nullopt; u.xcRadial=*n; }
        if (a==Axis::Hartree) { auto n=RungAbove(kAtomHartreeRadials, g.coulomb.nRadial); if (!n) return std::nullopt; u.coulomb.nRadial=*n; }
        if (a==Axis::Fit)     { if (g.fitRung+1>=BasisSet::Radial::CoulombFitLadder::kCoulombFitRungs) return std::nullopt; u.fitRung=g.fitRung+1; }
        return u;
    };
    std::vector<std::pair<Axis,bool>> axes{{Axis::XC, itsGrids.xcStated}};   // (axis, stated)
    if (itsGrids.coulomb.fit==H::CoulombFit::Multipole) axes.push_back({Axis::Hartree, itsGrids.hStated});
    if (dunlap && ladder)                               axes.push_back({Axis::Fit, false});

    for (int round=0; round<kAtomPilotRounds; round++)
    {
        rpt::Row row("rounds");
        const qchem::ChargeDensity::rDM_CD* d=itsDensity.get();
        auto key=[](const Grids& g){return std::to_string(g.xcRadial)+"/"+std::to_string(g.coulomb.nRadial)+"/"+std::to_string(g.fitRung);};
        std::map<std::string, std::vector<double>> probes;   // each grid sampled once at this D: every term's energy
        auto sample=[&](const H::rHamiltonian& h){ std::vector<double> E; for (const auto& te : h.MakeTermEnergies(d)) E.push_back(te.E); return E; };
        itsHam->RefreshForDensity(itsBasis, d);
        probes.emplace(key(itsGrids), sample(*itsHam));
        auto probe=[&](const Grids& g) -> const std::vector<double>&
        {
            if (auto it=probes.find(key(g)); it!=probes.end()) return it->second;
            // FENCED: a reference's own mesh and fit reports go to a nested run, so this run's `hamiltonian` keeps describing what it USED.
            struct Fence { explicit Fence(const std::string& n) {rpt::Begin(n);} ~Fence() {rpt::End();} } fence("atom pilot reference "+key(g));
            H::CoulombFitSlot<double> slot;   // a reference's own fit (a raised rung), never the run's
            const std::unique_ptr<H::rHamiltonian> h(MakeHamiltonianOn(g, slot));
            h->RefreshForDensity(itsBasis, d);
            return probes.emplace(key(g), sample(*h)).first->second;
        };
        auto distance=[&](const Grids& a, const Grids& b)
        {
            const auto& ea=probe(a); const auto& eb=probe(b);
            if (ea.size()!=eb.size()) throw std::logic_error("AtomCalculation::PilotGrids: two grids hand back different term lists");
            double worst=0.0;
            for (size_t t=0;t<ea.size();t++) worst=std::max(worst, std::fabs(ea[t]-eb[t]));
            return worst;
        };
        Grids next=itsGrids;
        bool raised=false;
        rpt::json steps=rpt::json::array();
        for (const auto& [axis, stated] : axes)
        {
            Grids cur=next;
            std::optional<double> dBelow;
            int raises=0;
            for (;;)
            {
                const auto up=above(cur, axis);
                if (!up)
                {
                    std::cout<<tag<<"WARNING: "<<name(axis)<<" "<<value(cur,axis)<<" is the ladder's top: unaudited, not settled"<<std::endl;
                    policy::Report("AtomGridLadderTop", name(axis)+" at its ladder's top "+std::to_string(value(cur,axis)), policy::Verdict::Warn,
                                   "no rung above it to measure against, so its error against eps "+std::to_string(target)+" is unknown");
                    break;
                }
                const double dd=distance(cur, *up), est=qcMesh::BeckeTailEstimate(dd, dBelow);
                const bool over=est>target;
                steps.push_back({{"axis",name(axis)}, {"rung",value(cur,axis)}, {"reference",value(*up,axis)}, {"delivered",dd},
                                 {"tailEstimate",est}, {"over",over}});
                std::cout<<tag<<"round "<<round<<"  "<<name(axis)<<" "<<value(cur,axis)<<" vs "<<value(*up,axis)<<": |ΔE| "
                         <<std::scientific<<std::setprecision(2)<<dd<<", tail "<<est<<(over ? " > " : " <= ")<<target<<std::defaultfloat<<std::endl;
                if (!over) break;
                if (stated)
                {
                    policy::Report("AtomStatedGridOverTarget", "the stated "+name(axis)+" "+std::to_string(value(cur,axis))+" misses eps",
                                   policy::Verdict::Warn, "a stated grid is audited, never raised (pin 31): tail estimate "+std::to_string(est)+
                                   " > target "+std::to_string(target)+"; the rung above, "+std::to_string(value(*up,axis))+", would be the next try");
                    break;
                }
                dBelow=dd; cur=*up; raises++; raised=true;
                if (round>0 && raises>=1) break;   // a confirming round raises once and re-converges
            }
            next=cur;
        }
        rpt::Set("round", round); rpt::Set("steps", steps);
        rpt::Set("grids", rpt::json{{"xcRadial",next.xcRadial}, {"hartreeRadial",next.coulomb.nRadial}, {"coulombFitRung",next.fitRung}});
        if (!raised) { std::cout<<tag<<"settled: every axis meets eps "<<target<<" at this D"<<std::endl; break; }
        if (round+1==kAtomPilotRounds)
        {
            std::cout<<tag<<"WARNING: the last round raised; its grids are NOT audited at the re-converged D"<<std::endl;
            policy::Report("AtomPilotRoundCap", "the atom pilot stops after "+std::to_string(kAtomPilotRounds)+" rounds", policy::Verdict::Warn,
                           "the last round raised a grid, and its re-converged density is not audited again");
        }
        // Raised: rebuild on the raised grids (a fresh fit slot when the rung moved) and re-converge ONCE from the carried D.
        if (next.fitRung!=itsGrids.fitRung) itsCoulombFit=std::make_unique<H::CoulombFitSlot<double>>();
        itsGrids=next;
        itsCarried=itsScf->GetWaveFunction()->GetChargeDensity();
        BuildGraph();
        if (itsObserver) itsScf->SetObserver(itsObserver);
        itsScf->Iterate(params);
        RebuildSampling();
        if (!itsScf->Converged()) { std::cout<<tag<<"WARNING: the re-convergence on the raised grids did not converge"<<std::endl; break; }
        if (round+1==kAtomPilotRounds) break;
    }
}

} // namespace qchem
