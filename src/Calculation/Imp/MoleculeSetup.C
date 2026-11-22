// File: Calculation/Imp/MoleculeSetup.C  The molecular run setup: every decision made ONCE, before a geometry point is built.
module;
#include <memory>
#include <cctype>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
#include <filesystem>
#include <vector>
module qchem.MoleculeSetup;

import qchem.PeriodicTable;   // thePeriodicTable (Z -> element symbol, for the PP lookup)
import qchem.Policy;          // policy::Refuse / Report: each policy files its row where it fires
import qchem.XCRhoPolicy;     // ReportXCRho: the XC ρ-axis decision
import qchem.ChargeDensity.AtomicDensity;   // HasAtomicDensity: the SAD seed's coverage
import qchem.SetupDecisions;  // the structure-neutral decisions (spin channels, seed, accelerator, Vxc admission), shared with SolidSetup
import qchem.BasisSet.Gaussian.Point.Factory; // RegisteredBasis, BasisSetData: which registered files are all-electron
import qchem.BasisSet.Gaussian.Point.BasisFiles; // BasisFileOfName: the basis name's file, for the coverage pre-flight
import qchem.BasisCoverage;   // RequireBasisCovers, BasisCovers, ElementsOf: a basis lacking an element is refused before anything is built
import qchem.Structure.MeshFactory;   // StructureMeshMaker<S>: each geometry point's mesh maker, typed where S is known (V-MESH-FACTORY)

namespace qchem
{

// THE MOLECULE'S VIEW OF THE ONE OPTIONS TREE: its per-kind defaults (each a `policy` row) and its refusals of what it does not consume.
RunOptions MoleculeOptions(const RunOptions& asked)
{
    RunOptions o=asked;
    const RunOptions def;
    auto refuse=[](const std::string& what, const std::string& why)
    { policy::Refuse("ForeignOption", what+" stated for a MOLECULE", why); };
    if (!(asked.atom==def.atom))   refuse("the `atom` block", "it holds the radial atom's options (its exponent pool, its pseudo-atom) -- drop it");
    if (!(asked.solid==def.solid))
        refuse("the `solid` block", "it holds periodic-only options a molecule cannot run -- drop it");
    if (asked.Nelec!=0)     refuse("Nelec", "a molecule counts its electrons from the structure and its charge -- drop it");
    if (asked.charge!=0)    refuse("charge", "a molecule's charge is its structure's -- state it on the structure");
    if (asked.ortho || asked.orthoTol!=0.0) refuse("ortho / orthoTol", "the molecular solver picks its own orthogonalisation (Auto)");
    if (!asked.spins.empty()) refuse("spins", "a per-site magnetic decoration is a periodic imposition's");
    if (!o.model) { o.model=Hamiltonian::Model::HF; policy::Report("DefaultModel", "model unstated: a molecule's default", policy::Verdict::Choose, "HF"); }
    if (!o.pp)     o.pp=std::string();                      // all-electron
    if (!o.xcMesh) o.xcMesh=kMolecularMeshDefault;          // climbed to eps from its entry rung
    return o;
}

// THE COULOMB FIT (policy CoulombFit; OpenWork "The Coulomb-fit default is a hidden policy", lane S 2026-11-21).  A density model's
// Hartree term on the Auxiliary route is a Dunlap fit on a Gaussian fit basis, and WHICH basis is a policy: before this row it was
// hard-wired in the orbital families (PG_Cart / PG_Spherical CreateCDFitBasisSet read A1_coul.bsd), against "policies live only in the
// policy layer".  Here the route's Auto resolves to Auxiliary, the fit file is named (kMolecularCoulombFitDefault when unstated), checked
// for existence and coverage, and MoleculeCalculation builds it into the run's CoulombFitSlot before the first Hamiltonian asks.
// The default is DGauss A1 (H-Te), today's numbers: it was built for DGauss's all-electron bases, and on a diffuse valence basis it is
// mHa low (DZVP-MOLOPT-GTH, LDA: OH- 1.9 mHa, H2O 4.4 mHa vs the delta route and PySCF's exact J).  The default per orbital-basis kind
// is the user's to rule on (step 3 of the row).
constexpr const char* kMolecularCoulombFitDefault="A1_coul";

static void ResolveCoulombFit(Hamiltonian::CoulombRoute& c, const Structure& point, bool densityModel)
{
    using CF=Hamiltonian::CoulombFit;
    const std::string stated=c.auxBasis;
    if (!densityModel)
    {
        if (!stated.empty())
            policy::Refuse("CoulombAuxBasisUnused", "coulomb.auxBasis \""+stated+"\" on a Hartree-Fock run",
                           "HF builds J from the exact 4-centre integrals and reads no Coulomb fit basis -- drop coulomb.auxBasis");
        return;
    }
    if (c.fit==CF::Delta || c.fit==CF::Multipole)
    {
        if (!stated.empty())
            policy::Refuse("CoulombAuxBasisUnused", "coulomb.auxBasis \""+stated+"\" on the "+(c.fit==CF::Delta ? "Delta" : "Multipole")+" Hartree route",
                           "that route puts rho on the XC mesh and has no fit basis -- drop coulomb.auxBasis, or state coulomb.fit Auxiliary");
        return;
    }
    const bool fitStated=c.fit==CF::Auxiliary;
    c.fit=CF::Auxiliary;
    if (stated.empty()) c.auxBasis=kMolecularCoulombFitDefault;
    const std::string file=c.auxBasis+".bsd";
    if (!std::filesystem::is_regular_file(BasisSet::Gaussian::BasisFile(file)))
    {
        std::string fits;
        for (const std::string& n : BasisSet::Gaussian::BasisNames())
            if (n.size()>5 && n.substr(n.size()-5)=="_coul") fits+=(fits.empty() ? "" : ", ")+n;
        policy::Refuse("CoulombAuxBasisUnknown", "coulomb.auxBasis \""+c.auxBasis+"\": no "+file+" in the basis data directory",
                       "state a Coulomb fit file by its stem: "+fits);
    }
    std::string missing;
    for (int Z : ElementsOf(point)) if (!BasisCovers(file, Z)) missing+=(missing.empty() ? "" : ", ")+thePeriodicTable().GetSymbol(Z);
    if (!missing.empty())
        policy::Refuse("CoulombAuxBasisCoverage", "Coulomb fit basis "+c.auxBasis+" has no block for "+missing,
                       std::string(stated.empty() ? "the default fit (A1, H-Te) does not cover this structure" : "the stated fit does not cover this structure")
                       +" -- state coulomb.fit Delta (rho on the XC mesh, every element), or a fit file that covers it");
    if (!stated.empty() && fitStated) return;   // fully stated: as asked, nothing decided
    policy::Report("CoulombFit", "the Hartree term: a Dunlap fit on "+c.auxBasis, policy::Verdict::Choose,
                   stated.empty()
                     ? std::string(fitStated ? "coulomb.fit Auxiliary stated without coulomb.auxBasis" : "no Hartree route stated")
                       +": the molecular default fit, DGauss A1 (H-Te, built for all-electron bases; mHa low on a diffuse valence basis such "
                       "as DZVP-MOLOPT-GTH) -- state coulomb.auxBasis (A1_coul, A2_coul) or coulomb.fit Delta to choose"
                     : "coulomb.auxBasis stated: its route is the Auxiliary (Dunlap) fit");
}

// THE DOORS (as SolidSetup's): Admit for user input, Resolve its throwing twin for code that states its own options.
Outcome<std::shared_ptr<const MoleculeSetup>, SetupRefusal> MoleculeSetup::Admit(const Molecule& st, const RunOptions& opts, const AcceleratorOptions& acc)
{
    return AdmitBy([&]{ return Resolve(st, opts, acc); });
}

// THE TYPED POINT BUILDER (V-MESH-FACTORY §8 step 2): made here, where the reference's concrete type S is known, so every geometry point
// is a copy of an S (never a Clone through the abstract face) and its mesh maker is StructureMeshMaker<S> over THAT point -- no cast.
// A pseudopotential run re-expresses the point as its VALENCE ions: each atom keeps its true species Z (so the PP lookup and the true
// geometry are intact) but ADDS the core charge Z-Zion (its own element, from the run's chosen set) to the net charge the site states,
// so GetNumElectrons() reports the Zion valence count minus the stated charge (OH- with GTH: 6 + 1 + 1), the count the EC + FittedVee
// charge constraint consume.
template <class S> static auto TypedPointBuilder(std::shared_ptr<const S> ref)
{
    return [ref](const std::vector<rvec3_t>* R, const Pseudopotential::PseudopotentialSet* pp)
    {
        auto point=std::make_shared<S>(*ref);
        if (R) point->MoveSites(*R);
        if (pp) point->AddNetCharges([pp](int Z){return double(Z - pp->Valence(thePeriodicTable().GetSymbol(Z)));});   // + the core Z - Zion: a stated ion stays one
        return MoleculeSetup::Point{point, std::make_shared<const StructureMeshMaker<S>>(point)};
    };
}

MoleculeSetup::Point MoleculeSetup::MakePoint() const {return itsPointBuilder(nullptr, pp ? &*pp : nullptr);}
MoleculeSetup::Point MoleculeSetup::MakePoint(const std::vector<rvec3_t>& R) const {return itsPointBuilder(&R, pp ? &*pp : nullptr);}

std::shared_ptr<const MoleculeSetup> MoleculeSetup::Resolve(const Molecule& st, const RunOptions& asked, const AcceleratorOptions& acc)
{
    auto s=std::make_shared<MoleculeSetup>();
    auto ref=std::make_shared<const Molecule>(st);   // deep copy, TRUE charges
    s->structure=ref;
    s->itsPointBuilder=TypedPointBuilder(ref);
    return Decide(std::move(s), asked, acc);
}

std::shared_ptr<const MoleculeSetup> MoleculeSetup::Decide(std::shared_ptr<MoleculeSetup> s, const RunOptions& asked, const AcceleratorOptions& acc)
{
    s->opts=MoleculeOptions(asked);   // the kind's defaults filled, what a molecule does not consume refused (V-FACADE-PARITY P0b-2)
    s->acc =acc;
    const RunOptions& opts=s->opts;   // from here on, the RESOLVED options
    // An open shell (2S>0) is unrestricted: it needs distinct up/down densities, so promote to polarized.
    s->opts.spin = ResolveSpinGroup(s->opts.spin, s->opts.multiplicity);   // THE rule, every structure (R6)
    // The pseudopotential front door builds a valence LDA Hamiltonian with the run's functional, LDA (Dirac + VWN5) or PADE (the
    // GTH-PADE PPs' own, 2026-11-12); any other model would be ignored, a SILENT substitution: refuse it, before anything is built.
    RouteMolecularGradient(*opts.model, opts.vxcFit);   // a GGA: the delta route samples grad rho (GGA stage 2); other fits refuse
    if (!opts.pp->empty() && *opts.model != Model::LDA && *opts.model != Model::PADE)
        policy::Refuse("PseudopotentialIsLDA", "a pseudopotential run with a model other than LDA or PADE",
                       "the molecular pseudopotential route builds a valence LDA Hamiltonian, so the model asked for would be "
                       "ignored -- state model LDA or PADE, or run all-electron (no pseudopotential)");
    if (!opts.pp->empty()) s->pp=Pseudopotential::PseudopotentialSet::Get(*opts.pp).With(opts.species);
    // An ALL-ELECTRON basis under a pseudopotential keeps core shells for electrons the pseudopotential has removed: the run is
    // computable (the core shells sit high and empty), but it is not the valence calculation the PP was generated for, and the
    // conditioning can suffer.  Known only for the registered files: an unregistered stem is not classified (no row).
    if (!opts.pp->empty())
    {
        using BasisSet::Gaussian::BasisSetData;
        BasisSetData data{};
        if (BasisSet::Gaussian::RegisteredBasis(opts.molecule.basis, data)
            && (data==BasisSetData::DZVP || data==BasisSetData::DZVP2 || data==BasisSetData::TZVP || data==BasisSetData::ORB
                || data==BasisSetData::ORB1))
            policy::Report("ValenceBasisUnderPP", "an all-electron basis ("+opts.molecule.basis+") under a pseudopotential", policy::Verdict::Warn,
                           "its core shells describe electrons the pseudopotential removed -- the run goes on; a valence basis "
                           "(sipp*, valence_lowq*) is the matched choice");
    }
    s->imposed=opts.imposeSymmetry;   // a molecule's imposition has no policy veto: as asked
    const std::shared_ptr<const Structure> point=s->MakePoint().structure;
    s->Ne=int(point->GetNumElectrons());
    RequireBasisStated(opts.molecule.basis, *point);   // no default basis: an unstated one is refused, naming the bases that cover this structure
    // Coverage, when the name is a basis file (an unknown name is the build's refusal, which lists every name there).
    if (!opts.molecule.basisFile.empty())
        RequireBasisCovers(opts.molecule.basisFile, opts.molecule.basis, *point);   // a stated file: it must cover every element (a missing one throws)
    else
    {
        std::string stem=opts.molecule.basis;
        for (char& c : stem) c=char(std::tolower((unsigned char)c));
        const auto names=BasisSet::Gaussian::BasisNames();
        if (std::find(names.begin(), names.end(), stem)!=names.end())
            RequireBasisCovers(BasisSet::Gaussian::BasisFileOfName(opts.molecule.basis), opts.molecule.basis, *point);
    }

    // The multiplicity 2S+1 -> (nUp, nDown): THE rule every structure shares (qchem.SetupDecisions), refusing an unrealisable one.
    // A molecule has no shared Fermi level to offer (nullopt).
    {
        const SpinChannels c=ResolveSpinChannels(opts.multiplicity, s->Ne, std::nullopt, "qchem::MoleculeCalculation");
        s->nUp=c.nUp; s->nDown=c.nDown;
    }

    // Seed: an explicit opts.seed wins; otherwise THE Auto rule both setups share (qchem.SeedPolicy: a model whose iteration-0 Fock reads a
    // density -- DFT, or any pseudopotential run -- seeds from the SAD superposition when the run's seed library has every element,
    // neutral here because the structure does not repeat, ionic on a lattice).  The library is the run's own atoms' (SeedKeyOf, as the
    // solid): a pseudopotential run reads the VALENCE densities at its q, an all-electron run the all-electron ones (before 2026-11-12 a
    // PP molecule read the all-electron densities, core included, and a species with no all-electron entry, Si, took the core guess).
    using qchem::ChargeDensity::SeedStrategy;
    const bool dft = qchem::Hamiltonian::IsDFT(*opts.model);
    {
        const RunSeed seed=ResolveRunSeed(opts.seed, *point, s->pp ? &*s->pp : nullptr, dft || !opts.pp->empty());   // no charge order stated
        s->opts.seed=seed.seed;
        s->seedKey=seed.key;
    }

    // THE XC MESH (pin 31; doc/GeometryOptimizationPlan.md §14.1): an axis left at 0 is unstated and takes the entry rung's; a mesh with
    // neither stated climbs from there at the converged D (MoleculeCalculation's pilot), only up.  A stated one is honoured and audited.
    s->xcMesh=*opts.xcMesh;
    const bool statedR=opts.xcMesh->nRadial>0, statedL=opts.xcMesh->angularDegree>0;
    s->meshStated=statedR || statedL;
    if (!statedR) s->xcMesh.nRadial      =kMolecularMeshEntryRadial;
    if (!statedL) s->xcMesh.angularDegree=kMolecularMeshEntryDegree;
    if (dft || !opts.pp->empty())
    {
        const std::string rung="("+std::to_string(s->xcMesh.nRadial)+", "+std::to_string(s->xcMesh.angularDegree)+")";
        if (!s->meshStated)
            policy::Report("MolecularXCMeshClimbs", "the Becke XC mesh climbs from "+rung+" to eps at the converged density", policy::Verdict::Choose,
                           "no XC mesh stated: the run sizes it from the accuracy target (pin 31), from the entry rung and only up");
        else if (!statedR || !statedL)
            policy::Report("MolecularXCMeshPartlyStated", "the Becke XC mesh "+rung+": the unstated axis takes the entry rung's", policy::Verdict::Choose,
                           std::string("the mesh states only its ")+(statedR ? "radial count" : "angular degree")+"; a stated mesh is audited, never raised");
    }

    // An HF run seeded from SAD reaches its first density matrix through ONE LDA (Xalpha) Fock (exact exchange needs a matrix, SAD has
    // none), built on THIS mesh: the entry rung when none is stated, the mesh a DFT user gets (lane A 2026-11-10: the old hidden
    // 12-direction rule split O2's pi_g by 13 mHa at the seed, and a stated Shared fill converged to the broken state).
    if (!dft && opts.pp->empty() && s->opts.seed==SeedStrategy::SAD)
        policy::Report("HFSeedFockMesh", "the SAD seed reaches HF through one LDA Fock on the XC mesh ("+std::to_string(s->xcMesh.nRadial)+", "
                       +std::to_string(s->xcMesh.angularDegree)+")", policy::Verdict::Choose,
                       std::string(s->meshStated ? "the stated mesh" : "the molecular entry rung, no mesh stated")
                       +": a seed Fock coarser than the run's mesh can break the molecule's symmetry before the first iteration");

    // The Vxc product space's admission (pin 30) and the Hartree route's, asked of the Hamiltonian's rules HERE, before the basis
    // and the Hamiltonian are built: the Factory throws the same reason, but only once the work is done.
    if (dft || !opts.pp->empty())
    {
        RequireVxcAdmitted(*point, opts.vxcFit, s->xcMesh, "molecule");
        if (const auto a=Hamiltonian::AdmitCoulomb(*point, opts.coulomb, s->xcMesh, opts.vxcFit); !a)
            policy::Refuse("AdmitCoulomb", "the stated Hartree route on this molecule", a.Error());
        ReportXCRho(opts.vxcFit, opts.xcRho);   // the ρ axis (user ruling 2026-11-09): warns on the Gaussian x δ cell
    }
    ResolveCoulombFit(s->opts.coulomb, *point, dft || !opts.pp->empty());   // policy CoulombFit: the Hartree fit named here, never in a library

    AdmitRecipe(opts.recipe, /*kT: not read for a molecule*/0.0, /*molecule*/true);   // Auto | Insulator only (R4.1 ruling (c))
    s->accelerator = MolecularAcceleratorType(opts.accelerator, "qchem::MoleculeCalculation");   // Auto = DIIS; Null refused
    if (dft || !opts.pp->empty()) RefuseGradientConsumers(opts.vxcFit, opts.xcRho, NameOf(s->accelerator));   // GDM needs H = dE/dD (pin 32)
    const double eMax=DIISEngageGate(acc.diis.engageBelowFD.value_or(0.0), point->GetNuclearCharge(), dft || !opts.pp->empty());   // DIIS from the start on a density Fock
    s->accParams = ResolveAcceleratorParams(acc, MolecularAcceleratorDefaults(eMax));   // the gate lands in accParams.diis.EngageBelowFD
    return s;
}

std::vector<Hamiltonian::PPFunctionalCheck> MoleculeSetup::PPFunctionalChecks() const
{
    if (!pp || !Hamiltonian::IsDFT(*opts.model)) return {};
    return Hamiltonian::PPFunctionalChecks(*pp, *structure, Hamiltonian::ModelToXC(*opts.model, opts.xalpha));
}


} // namespace qchem
