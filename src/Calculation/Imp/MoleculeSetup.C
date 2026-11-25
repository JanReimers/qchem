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
#include <sstream>   // the MeshSymmetry row's numbers
#include <iomanip>
#include <cmath>
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
import qchem.Symmetry.Molecule.PointGroup;     // DetectPointGroup: the point group the mesh's policy names (policy MeshSymmetry)
import qchem.Symmetry.Molecule.AbelianGroup;   // BuildAbelianGroup: its D2h-type subgroup, the ops the mesh is made invariant under

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

// THE COULOMB FIT (policy CoulombFit; OpenWork "The Coulomb-fit default is a hidden policy", lane S 2026-11-21/22).  A density model's
// Hartree route and, on the Auxiliary route, WHICH Gaussian fit basis are policies: before this row A1 was hard-wired in the orbital
// families (PG_Cart / PG_Spherical CreateCDFitBasisSet read A1_coul.bsd), against "policies live only in the policy layer".
// USER RULING 2026-11-22 ("Yes, δ default"):
//  - an unstated route is DELTA (rho on the XC mesh, radial Poisson) for every Gaussian orbital basis: the one route that meets eps.
//    Measured (lane-s-CoulombFit-step3.md): the DGauss fits miss eps on every molecule tried, DGauss's own all-electron dzvp included --
//    A1 0.4-5 mHa, A2 0.08-1.2 mHa on Cartesian orbitals, up to 0.1 Ha on spherical ones (H2S) -- while delta is 0.4 uHa from PySCF's
//    exact J (OH-).  Its price is the mesh climb its mesh-borne Hartree needs: 2.2x (TeF6) to 8.6x (SF6) an A fit's wall.
//  - Auxiliary stated without a file takes A2 (kMolecularCoulombFitDefault), the best Gaussian fit measured; where A2 does not cover the
//    structure (it stops at Zn and has no K, Ca) A1 (H-Te) instead, said so.
//  - any stated file runs, on Cartesian or spherical orbitals; MoleculeCalculation's CoulombFitResidual check WARNs when it misses eps.
// MoleculeCalculation builds a named fit into the run's CoulombFitSlot before the first Hamiltonian asks.
constexpr const char* kMolecularCoulombFitDefault ="A2_coul";
constexpr const char* kMolecularCoulombFitFallback="A1_coul";

static void ResolveCoulombFit(Hamiltonian::CoulombRoute& c, const Structure& point, bool densityModel)
{
    using CF=Hamiltonian::CoulombFit;
    const std::string stated=c.auxBasis;
    if (!densityModel)
    {
        if (!stated.empty())
            policy::Refuse("CoulombAuxBasisUnused", "coulomb.auxBasis \""+stated+"\" on a Hartree-Fock run",
                           "HF builds J from the exact 4-centre integrals and reads no Coulomb fit basis -- drop coulomb.auxBasis");
        if (c.metric!=Hamiltonian::CoulombMetric::Auto || c.constraint!=Hamiltonian::ChargeConstraint::Auto)
            policy::Refuse("CoulombFitAxesUnused", "coulomb.metric / coulomb.constraint stated on a Hartree-Fock run",
                           "HF builds J from the exact 4-centre integrals and fits nothing -- drop coulomb.metric and coulomb.constraint");
        return;
    }
    if (c.fit==CF::Delta || c.fit==CF::Multipole)
    {
        if (!stated.empty())
            policy::Refuse("CoulombAuxBasisUnused", "coulomb.auxBasis \""+stated+"\" on the "+(c.fit==CF::Delta ? "Delta" : "Multipole")+" Hartree route",
                           "that route puts rho on the XC mesh and has no fit basis -- drop coulomb.auxBasis, or state coulomb.fit Auxiliary");
        return;
    }
    // A stated fit axis (V-FIT-AXES) implies the Gaussian route, as a stated fit file does: only that fit has the axes to choose.
    const bool axesStated=c.metric!=Hamiltonian::CoulombMetric::Auto || c.constraint!=Hamiltonian::ChargeConstraint::Auto;
    if (c.fit==CF::Auto && stated.empty() && !axesStated)
    {
        c.fit=CF::Delta;
        policy::Report("CoulombFit", "the Hartree term: delta (rho on the XC mesh, radial Poisson)", policy::Verdict::Choose,
                       "no Hartree route stated: the molecular default (user ruling 2026-11-22), the one route that meets eps -- the DGauss "
                       "fits are mHa low on every molecule measured; state coulomb.fit Auxiliary (A2) or coulomb.auxBasis for the faster fit");
        return;
    }
    const bool fitStated=c.fit==CF::Auxiliary;
    c.fit=CF::Auxiliary;
    // THE FIT AXES (V-FIT-AXES): unstated, Dunlap's fit -- the Coulomb metric, charge-constrained (AdmitCoulomb refused what is not built).
    if (c.metric==Hamiltonian::CoulombMetric::Auto)        c.metric=Hamiltonian::CoulombMetric::Coulomb;
    if (c.constraint==Hamiltonian::ChargeConstraint::Auto) c.constraint=Hamiltonian::ChargeConstraint::Charge;
    const std::string axes=std::string(" (Coulomb metric, ")+(c.constraint==Hamiltonian::ChargeConstraint::Charge ? "charge-constrained" : "unconstrained")
                           +(axesStated ? ", as stated)" : ": the Dunlap default; state coulomb.constraint None for the plain fit)");
    auto uncovered=[&](const std::string& file)
    {
        std::string missing;
        for (int Z : ElementsOf(point)) if (!BasisCovers(file, Z)) missing+=(missing.empty() ? "" : ", ")+thePeriodicTable().GetSymbol(Z);
        return missing;
    };
    std::string fallbackWhy;
    if (stated.empty())
    {
        c.auxBasis=kMolecularCoulombFitDefault;
        if (const std::string m=uncovered(c.auxBasis+".bsd"); !m.empty())
        {
            fallbackWhy="A2 has no block for "+m+", so A1 (H-Te) instead";
            c.auxBasis=kMolecularCoulombFitFallback;
        }
    }
    const std::string file=c.auxBasis+".bsd";
    if (!std::filesystem::is_regular_file(BasisSet::Gaussian::BasisFile(file)))
    {
        std::string fits;
        for (const std::string& n : BasisSet::Gaussian::BasisNames())
            if (n.size()>5 && n.substr(n.size()-5)=="_coul") fits+=(fits.empty() ? "" : ", ")+n;
        policy::Refuse("CoulombAuxBasisUnknown", "coulomb.auxBasis \""+c.auxBasis+"\": no "+file+" in the basis data directory",
                       "state a Coulomb fit file by its stem: "+fits);
    }
    if (const std::string missing=uncovered(file); !missing.empty())
        policy::Refuse("CoulombAuxBasisCoverage", "Coulomb fit basis "+c.auxBasis+" has no block for "+missing,
                       std::string(stated.empty() ? "neither A2 nor A1 (H-Te) covers this structure" : "the stated fit does not cover this structure")
                       +" -- state coulomb.fit Delta (rho on the XC mesh, every element), or a fit file that covers it");
    if (!stated.empty() && fitStated) return;   // fit and file stated: as asked, nothing decided (unstated axes are the long-standing Dunlap default)
    const std::string route = fitStated ? "" : (stated.empty() ? "a stated fit axis implies the Auxiliary route; " : "coulomb.auxBasis stated: its route is the Auxiliary fit; ");
    policy::Report("CoulombFit", "the Hartree term: a Gaussian density fit on "+c.auxBasis+axes, policy::Verdict::Choose,
                   route+(stated.empty()
                     ? "no fit file stated: the fastest-good fit, DGauss A2 (user ruling 2026-11-22)"
                       +(fallbackWhy.empty() ? std::string() : " -- "+fallbackWhy)+"; its residual is checked against eps (CoulombFitResidual)"
                     : "the stated file "+c.auxBasis));
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

MoleculeSetup::Point MoleculeSetup::MakePoint() const
{
    Point p=itsPointBuilder(nullptr, pp ? &*pp : nullptr);
    p.meshSymmetry=MeshSymmetryAt(*p.structure, /*report*/false);
    return p;
}
MoleculeSetup::Point MoleculeSetup::MakePoint(const std::vector<rvec3_t>& R) const
{
    Point p=itsPointBuilder(&R, pp ? &*pp : nullptr);
    p.meshSymmetry=MeshSymmetryAt(*p.structure, /*report*/false);
    return p;
}

namespace {
//! A molecule's point-group mesh source: the ops of its D2h-type subgroup H, {W | 0} about the centroid (the SymOp point-group
//! convention; the maker supplies the centre), copied from the SALC group's ops -- the adapter V-SYM-ONE-SYMOP will retire.
class MolecularPointGroupSource : public qcMesh::SymmetrySource
{
public:
    MolecularPointGroupSource(std::vector<Symmetry::SymOp> ops, bool averages) : itsOps(std::move(ops)), itsAverages(averages) {}
    std::vector<Symmetry::SymOp> MeshOps() const override {return itsOps;}
    bool AveragesDensityOverOrbits() const override {return itsAverages;}
private:
    std::vector<Symmetry::SymOp> itsOps;
    bool                         itsAverages;
};
} //anon

// POLICY MeshSymmetry (OOD-SOLID V-FOLD-NEUTRAL; user ruling 2026-11-23: D1 the abelian SALC subgroup H now, D2 free runs too, D3 a C1
// molecule stays in the lab frame).  The mesh is made invariant under H -- the one shared angular set rotated into H's frame, which
// removes the deck's orientation from the answer -- and folded when the run imposes symmetry.  Three cases keep the lab frame, each a
// reported choice: an angular rule that cannot carry H (Lebedev below degree 9, GaussLegendre), and a geometry that carries H only to
// the detection tolerance (the oriented mesh needs kMolecularMeshGeometryTol: its members' grids sit on their own nuclei).
std::shared_ptr<const qcMesh::SymmetrySource> MoleculeSetup::MeshSymmetryAt(const Structure& point, bool report) const
{
    auto none=std::make_shared<const qcMesh::NoSymmetry>();
    const bool dft=Hamiltonian::IsDFT(*opts.model);
    if (!dft && opts.pp->empty()) return none;   // no XC mesh: HF all-electron (its SAD seed Fock's mesh stays plain)
    namespace SM=Symmetry::Molecule;
    std::vector<SM::SymPoint> pts;
    point.ForEachSite([&pts](int Z, const rvec3_t& r, bool){ pts.push_back({Z, r}); });
    const SM::PointGroup pg=SM::DetectPointGroup(pts, opts.symmetryTol);
    if (pg.symbol=="C1") return none;            // D3: no group, no frame to orient to
    const SM::AbelianGroup H=SM::BuildAbelianGroup(pts, opts.symmetryTol);
    const std::string group=pg.symbol+(pg.abelian==pg.symbol ? std::string() : " (its D2h-type subgroup "+pg.abelian+")");
    auto row=[&](const std::string& what, const std::string& why)
    {
        if (report) policy::Report("MeshSymmetry", what, policy::Verdict::Choose, why);
        else if (group!=itsMeshGroup)   // a later point: only a change of group
            policy::Report("MeshSymmetry", what+" (this geometry point; the reference's group was "+(itsMeshGroup.empty() ? "C1" : itsMeshGroup)+")",
                           policy::Verdict::Choose, why);
    };
    if (H.ops.size()<=1) return none;            // a group with no D2h-type element beyond E (C3, ...): nothing to orient to
    // The angular rule must carry the frame's reflections (MolecularPointGroupMesh.EveryLebedevRuleFromDegree9CarriesTheFrame).
    if (xcMesh.angular!=qcMesh::AngularKind::Lebedev || xcMesh.angularDegree<9)
    {
        row("point group "+group+" detected; the XC mesh stays in the lab frame",
            "its angular rule (degree "+std::to_string(xcMesh.angularDegree)+") cannot carry the group's frame reflections: only a Lebedev rule of "
            "degree 9 or more does -- the answer then depends on how the deck orients the molecule, at the mesh's orientation error");
        return none;
    }
    // A ROTATED angular rule (angRot, the cubic-site calibration that steers Lebedev's special orbits off a crystal's bonds) is a stated
    // orientation of its own: the group's frame would undo it, so it is kept, and the mesh with it in the lab frame.
    if (xcMesh.angRot!=0.0)
    {
        std::ostringstream a; a << xcMesh.angRot;
        row("point group "+group+" detected; the XC mesh stays in the lab frame",
            "its angular rule is stated ROTATED (angRot "+a.str()+" rad), an orientation of its own that the group's frame would undo -- the "
            "answer then depends on how the deck orients the molecule, at the mesh's orientation error.  To fix this, please state angRot 0");
        return none;
    }
    // The geometry must carry H to the mesh's precision, not just to the detection tolerance (measured after the maker's polish).
    std::vector<Symmetry::SymOp> ops;
    for (const auto& op : H.ops) ops.push_back({op.Matrix(), rvec3_t(0,0,0)});
    const double d=PointGroupGeometryError(point, ops);
    if (d>kMolecularMeshGeometryTol)
    {
        std::ostringstream o; o << std::scientific << std::setprecision(1) << d;
        row("point group "+group+" detected; the XC mesh stays in the lab frame",
            "the geometry carries the group only to "+o.str()+" bohr, and an invariant mesh needs "+[]{std::ostringstream t; t<<std::scientific
            <<std::setprecision(0)<<kMolecularMeshGeometryTol; return t.str();}()+" -- the answer then depends on how the deck orients the "
            "molecule, at the mesh's orientation error.  To fix this, please state the coordinates symmetric to more digits");
        return none;
    }
    // THE FOLD (OpenWork §3, the per-orbit consumer, lane H 2026-11-24): engaged on an IMPOSED DFT run whose Vxc fit is δ.  There the basis is
    // SALC-blocked into 1-D real irreps, so every XC integrand (ρ, v χ_iχ_j, |∇ρ|², ∇ρ·∇(χ_iχ_j)) is invariant pointwise and the XC term
    // samples ONE point per orbit at the orbit's summed weight -- exact, not an approximation (FoldAudit, to 1e-12), with the force through
    // the A1 projection (HamiltonianForces.TheFoldedXCForceIsTheUnfoldedOne).  A free run's ρ is not invariant; a Gaussian (A1) or multipole
    // Vxc fit's functions are not invariant (their projection needs the orbit-summed table, a follow-up), and HF has no XC term to fold.
    const bool deltaFit=opts.vxcFit==Hamiltonian::VxcFit::Auto || opts.vxcFit==Hamiltonian::VxcFit::Delta;
    const bool fold=imposed && dft && deltaFit;
    const std::string frame="the one angular set is rotated into the group's frame (reported as nonlinearMesh.frame), so the deck's orientation "
        "does not reach the answer"+std::string(pg.abelian==pg.symbol ? "" : "; invariance under the full group "+pg.symbol+" is the follow-up "
        "(OOD-SOLID V-FOLD-NEUTRAL)");
    std::string what="point group "+group+" detected; the XC mesh is made invariant under "+pg.abelian+" ("+std::to_string(ops.size())+" ops)";
    std::string why=frame;
    if (fold)
    {
        what+=", FOLDED: the XC term samples one point per orbit (counted as perOrbit in the XC mesh's report)";
        why+="; on an imposed run every XC integrand is invariant under the group (a SALC basis of 1-D real irreps), so one point per orbit at "
             "the orbit's summed weight gives the same energy, matrix and force as the whole mesh, at about 1/|group| of the collocation";
    }
    else if (!imposed) what+=", not folded (a free run: its density is not invariant under the group)";
    else
    {
        what+=", NOT folded on this imposed run";
        why+=!dft ? std::string("; HF has no XC term to sample per orbit")
                  : std::string("; its Vxc fit is not the delta fit, and a Gaussian or multipole fit function is not invariant under the group, so "
                                "one point per orbit would mis-project it (the orbit-summed fit table is a follow-up, OpenWork §3) -- state the delta "
                                "Vxc fit (the default) for the per-orbit XC");
    }
    row(what, why);
    return std::make_shared<const MolecularPointGroupSource>(std::move(ops), /*averages*/fold);
}

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
    // The reference point's STRUCTURE only (its mesh symmetry is decided below, once the mesh recipe is resolved).
    const std::shared_ptr<const Structure> point=s->itsPointBuilder(nullptr, s->pp ? &*s->pp : nullptr).structure;
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

    // THE MESHES' SYMMETRY SOURCE (policy MeshSymmetry): decided per geometry point by MakePoint; the reference's files its row here, and
    // its group is kept so a later point reports only a CHANGE of group.
    {
        std::vector<Symmetry::Molecule::SymPoint> pts;
        point->ForEachSite([&pts](int Z, const rvec3_t& r, bool){ pts.push_back({Z, r}); });
        const Symmetry::Molecule::PointGroup pg=Symmetry::Molecule::DetectPointGroup(pts, opts.symmetryTol);
        (void)s->MeshSymmetryAt(*point, /*report*/true);
        s->itsMeshGroup=pg.symbol=="C1" ? std::string() : pg.symbol+(pg.abelian==pg.symbol ? std::string() : " (its D2h-type subgroup "+pg.abelian+")");
    }

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
