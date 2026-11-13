// File: Calculation/Imp/MoleculeSetup.C  The molecular run setup: every decision made ONCE, before a geometry point is built.
module;
#include <memory>
#include <cctype>
#include <algorithm>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>
module qchem.MoleculeSetup;

import qchem.PeriodicTable;   // thePeriodicTable (Z -> element symbol, for the PP lookup)
import qchem.Policy;          // policy::Refuse / Report: each policy files its row where it fires
import qchem.XCRhoPolicy;     // ReportXCRho: the XC ρ-axis decision
import qchem.ChargeDensity.AtomicDensity;   // HasAtomicDensity: the SAD seed's coverage
import qchem.SeedPolicy;      // ResolveSeed: the ONE Auto seed rule (shared with SolidSetup)
import qchem.BasisSet.Gaussian.Point.Factory; // RegisteredBasis, BasisSetData: which registered files are all-electron
import qchem.BasisSet.Gaussian.Point.BasisFiles; // BasisFileOfName: the basis name's file, for the coverage pre-flight
import qchem.BasisCoverage;   // RequireBasisCovers: a basis lacking an element is refused before anything is built

namespace qchem
{

std::shared_ptr<const MoleculeSetup> MoleculeSetup::Resolve(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc)
{
    auto s=std::make_shared<MoleculeSetup>();
    s->opts=opts;
    s->acc =acc;
    // An open shell (2S>0) is unrestricted: it needs distinct up/down densities, so promote to polarized.
    if (s->opts.multiplicity > 1) s->opts.spin = SpinGroup::Polarized;
    // The pseudopotential front door builds a valence LDA Hamiltonian with the run's functional, LDA (Dirac + VWN5) or PADE (the
    // GTH-PADE PPs' own, 2026-11-12); any other model would be ignored, a SILENT substitution: refuse it, before anything is built.
    if (!opts.pp.empty() && opts.model != Model::LDA && opts.model != Model::PADE)
        policy::Refuse("PseudopotentialIsLDA", "a pseudopotential run with a model other than LDA or PADE",
                       "the molecular pseudopotential route builds a valence LDA Hamiltonian, so the model asked for would be "
                       "ignored -- state model LDA or PADE, or run all-electron (no pseudopotential)");
    if (!opts.pp.empty()) s->pp=Pseudopotential::PseudopotentialSet::Get(opts.pp).With(opts.ppOverride);
    // An ALL-ELECTRON basis under a pseudopotential keeps core shells for electrons the pseudopotential has removed: the run is
    // computable (the core shells sit high and empty), but it is not the valence calculation the PP was generated for, and the
    // conditioning can suffer.  Known only for the registered files: an unregistered stem is not classified (no row).
    if (!opts.pp.empty())
    {
        using BasisSet::Gaussian::BasisSetData;
        BasisSetData data{};
        if (BasisSet::Gaussian::RegisteredBasis(opts.basis, data)
            && (data==BasisSetData::DZVP || data==BasisSetData::DZVP2 || data==BasisSetData::TZVP || data==BasisSetData::ORB
                || data==BasisSetData::ORB1))
            policy::Report("ValenceBasisUnderPP", "an all-electron basis ("+opts.basis+") under a pseudopotential", policy::Verdict::Warn,
                           "its core shells describe electrons the pseudopotential removed -- the run goes on; a valence basis "
                           "(sipp*, valence_lowq*) is the matched choice");
    }
    s->reference=st.Clone();   // polymorphic deep copy: a UnitCell stays periodic, not sliced to a Molecule
    const std::shared_ptr<const Structure> point=s->PointStructure(*s->reference);
    s->Ne=int(point->GetNumElectrons());
    RequireBasisStated(opts.basis, *point);   // no default basis: an unstated one is refused, naming the bases that cover this structure
    // Coverage, when the name is a basis file (an unknown name is the build's refusal, which lists every name there).
    if (!opts.basisFile.empty())
        RequireBasisCovers(opts.basisFile, opts.basis, *point);   // a stated file: it must cover every element (a missing one throws)
    else
    {
        std::string stem=opts.basis;
        for (char& c : stem) c=char(std::tolower((unsigned char)c));
        const auto names=BasisSet::Gaussian::BasisNames();
        if (std::find(names.begin(), names.end(), stem)!=names.end())
            RequireBasisCovers(BasisSet::Gaussian::BasisFileOfName(opts.basis), opts.basis, *point);
    }

    // The multiplicity 2S+1 -> (nUp, nDown) must be realisable: 0 <= 2S <= Ne and Ne-2S even (fail loud, not a silent miscount).
    if (opts.multiplicity > 0)
    {
        const int twoS = opts.multiplicity - 1;
        if (twoS > s->Ne || ((s->Ne - twoS) % 2) != 0)
            throw std::runtime_error("qchem::MoleculeCalculation: multiplicity " + std::to_string(opts.multiplicity) +
                " (2S=" + std::to_string(twoS) + ") is incompatible with " + std::to_string(s->Ne) +
                " electrons -- need 0 <= 2S <= Ne and Ne-2S even.");
    }

    // Seed: an explicit opts.seed wins; otherwise THE Auto rule both setups share (qchem.SeedPolicy: a model whose iteration-0 Fock reads a
    // density -- DFT, or any pseudopotential run -- seeds from the SAD superposition when the all-electron library has every element,
    // neutral here because the structure does not repeat, ionic on a lattice).  (Measured 2026-11-07: H2O LDA + GTH on dzvp converges in 9
    // iterations from SAD, 26 from the core guess, although the SAD densities are all-electron; Si has no LDA entry, so its PP runs take
    // the core guess, with a WARN since 2026-11-12.)
    using qchem::ChargeDensity::SeedStrategy;
    const bool dft = qchem::Hamiltonian::IsDFT(opts.model);
    s->seed = ResolveSeed(opts.seed, *point, {.densityFock = dft || !opts.pp.empty(), .periodic = point->HasTranslationSymmetry(),
                                             .chargesStated = false, .covered = [](int Z){ return ChargeDensity::HasAtomicDensity(Z); }});

    // THE XC MESH (pin 31; doc/GeometryOptimizationPlan.md §14.1): an axis left at 0 is unstated and takes the entry rung's; a mesh with
    // neither stated climbs from there at the converged D (MoleculeCalculation's pilot), only up.  A stated one is honoured and audited.
    s->mesh=opts.mesh;
    const bool statedR=opts.mesh.nRadial>0, statedL=opts.mesh.angularDegree>0;
    s->meshStated=statedR || statedL;
    if (!statedR) s->mesh.nRadial      =kMolecularMeshEntryRadial;
    if (!statedL) s->mesh.angularDegree=kMolecularMeshEntryDegree;
    if (dft || !opts.pp.empty())
    {
        const std::string rung="("+std::to_string(s->mesh.nRadial)+", "+std::to_string(s->mesh.angularDegree)+")";
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
    if (!dft && opts.pp.empty() && s->seed==SeedStrategy::SAD)
        policy::Report("HFSeedFockMesh", "the SAD seed reaches HF through one LDA Fock on the XC mesh ("+std::to_string(s->mesh.nRadial)+", "
                       +std::to_string(s->mesh.angularDegree)+")", policy::Verdict::Choose,
                       std::string(s->meshStated ? "the stated mesh" : "the molecular entry rung, no mesh stated")
                       +": a seed Fock coarser than the run's mesh can break the molecule's symmetry before the first iteration");

    // The Vxc product space's admission (pin 30) and the Hartree route's, asked of the Hamiltonian's rules HERE, before the basis
    // and the Hamiltonian are built: the Factory throws the same reason, but only once the work is done.
    if (dft || !opts.pp.empty())
    {
        if (const auto a=Hamiltonian::AdmitVxc(*point, opts.vxcFit, s->mesh); !a)
            policy::Refuse("AdmitVxc", "the stated Vxc fit and XC mesh on this molecule", a.Error());
        if (const auto a=Hamiltonian::AdmitCoulomb(*point, opts.coulomb, s->mesh, opts.vxcFit); !a)
            policy::Refuse("AdmitCoulomb", "the stated Hartree route on this molecule", a.Error());
        ReportXCRho(opts.vxcFit, opts.xcRho);   // the ρ axis (user ruling 2026-11-09): warns on the Gaussian x δ cell
    }

    using AType = qchem::SCFAccelerators::Type;
    s->accelerator = acc.type=="DIIS"   ? AType::DIIS
                   : acc.type=="GDM"    ? AType::GDM
                   : acc.type=="Ladder" ? AType::Ladder
                   : throw std::runtime_error("qchem::MoleculeCalculation: unknown accelerator type \"" + acc.type +
                                              "\" (expected DIIS | GDM | Ladder)");
    if (dft || !opts.pp.empty()) RefuseGradientConsumers(opts.vxcFit, opts.xcRho, acc.type);   // GDM needs H = dE/dD (pin 32)
    // DFT (and the LSDA pseudopotential) need DIIS engaged from iteration 0: the default Z-scaled EMax gate keeps DIIS off and
    // the SCF limit-cycles to a non-converged snapshot (the "M_Sym layout UB" mechanism, see M_DFT).  So default them to a high
    // EMax ("DIIS from start") unless the caller pinned one.
    const bool   dftLike = dft || !opts.pp.empty();
    const int    Z       = point->GetNuclearCharge();
    s->eMax = acc.eMax > 0.0 ? acc.eMax : (dftLike ? 100.0 : Z*Z*0.08/32);
    return s;
}

// Re-express a structure as its VALENCE ions for a pseudopotential run: each atom keeps its true species Z (so the PP lookup and
// the true geometry are intact) but carries net charge Z-Zion(its own element, from the run's chosen set), so GetNumElectrons()
// reports the total Zion valence count the EC + FittedVee charge constraint consume.  Clone-and-mutate (not rebuild-as-Molecule)
// so the concrete geometry is preserved: a PP run on a periodic UnitCell stays a UnitCell (Ewald ion-ion, uniform mesh).
std::vector<Hamiltonian::PPFunctionalCheck> MoleculeSetup::PPFunctionalChecks() const
{
    if (!pp || !Hamiltonian::IsDFT(opts.model)) return {};
    return Hamiltonian::PPFunctionalChecks(*pp, *reference, Hamiltonian::ModelToXC(opts.model, opts.xalpha));
}

std::shared_ptr<const Structure> MoleculeSetup::PointStructure(const Structure& st) const
{
    auto point = st.Clone();
    if (pp) point->SetNetCharges([this](int Z){return double(Z - pp->Valence(thePeriodicTable().GetSymbol(Z)));});   // net charge = Z - Zion
    return point;
}

nlohmann::json ToAcceleratorJson(const AcceleratorOptions& acc, double eMax)
{
    return {{"NProj",acc.nProj},{"EMax",eMax},{"EMin",acc.eMin},{"SVTol",acc.svTol},{"SVTolRel",acc.svTolRel},{"Trust",acc.trust},
            {"ethresh",acc.ladderEThresh},{"stall",acc.ladderStall},{"floor",acc.ladderFloor},{"switchat",acc.ladderSwitchAt},
            {"signal",acc.ladderSignal==SCFAccelerators::ScheduleSignal::EnergyChange ? "EnergyChange" : "Error"}};
}

} // namespace qchem
