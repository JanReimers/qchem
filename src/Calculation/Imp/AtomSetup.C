// File: Calculation/Imp/AtomSetup.C  The radial atom's run setup: every decision made ONCE, before its graph is built.
module;
#include <cmath>
#include <memory>
#include <string>
#include <nlohmann/json.hpp>
module qchem.AtomSetup;

import qchem.Policy;            // policy::Refuse / Report: each policy files its row where it fires
import qchem.SetupDecisions;    // AcceleratorTypeByName, DIISEngageGate, ResolveAcceleratorParams: the decisions every setup shares
import qchem.XCRhoPolicy;       // RefuseGradientConsumers (pin 32)
import qchem.ChargeDensity.Seed;   // SeedStrategy
import qchem.Types;                // rvec3_t (the atom at the origin)

namespace qchem
{

// THE ATOM'S VIEW OF THE ONE OPTIONS TREE (V-FACADE-PARITY P0b-2): its defaults and its refusals of what it does not consume.
RunOptions AtomOptions(const RunOptions& asked)
{
    RunOptions o=asked;
    const RunOptions def;
    auto refuse=[](const std::string& what, const std::string& why)
    { policy::Refuse("ForeignOption", what+" stated for an ATOM", why); };
    if (!(asked.molecule==def.molecule)) refuse("the `molecule` block", "an atom's basis is its radial exponent pool (the `atom` block)");
    if (!(asked.solid==def.solid))       refuse("the `solid` block", "it holds periodic-only options an atom cannot run");
    if (asked.Nelec!=0 || asked.charge!=0 || asked.multiplicity!=0)
        refuse("Nelec / charge / multiplicity", "an atom's Z and charge are its structure's and its spin is the configuration's");
    if (asked.pp || !asked.species.empty())
        refuse("pp / species", "an atom states its pseudo-atom in the `atom` block (pseudopotential, valence)");
    if (asked.imposeSymmetry || asked.symmetryTol!=def.symmetryTol || !asked.spins.empty())
        refuse("imposeSymmetry / symmetryTol / spins", "an atom's symmetry is spherical by construction");
    if (asked.orthoTol!=def.orthoTol) refuse("orthoTol", "the atom's solver takes no truncation floor");
    if (asked.label!=def.label) refuse("label", "an atom names itself by its element in its own output");
    if (!o.model)  o.model=Hamiltonian::Model::HF;
    if (!o.xcMesh) o.xcMesh=kAtomMeshDefault;
    if (!o.ortho)  o.ortho=qchem::Auto;
    return o;
}

Outcome<std::shared_ptr<const AtomSetup>, SetupRefusal> AtomSetup::Admit(const Structure& atom, const RunOptions& opts, const AcceleratorOptions& acc)
{
    return AdmitBy([&]{ return Resolve(atom, opts, acc); });
}

std::shared_ptr<const AtomSetup> AtomSetup::Resolve(int Z, int charge, const RunOptions& opts, const AcceleratorOptions& acc)
{
    return Resolve(Atom(Z, double(charge)), opts, acc);
}

std::shared_ptr<const AtomSetup> AtomSetup::Resolve(const Structure& atom, const RunOptions& asked, const AcceleratorOptions& acc)
{
    if (atom.GetNumAtoms()!=1)
        policy::Refuse("AtomStructure", "a structure of "+std::to_string(atom.GetNumAtoms())+" atoms for the radial atom",
                       "the radial atom solves ONE atom at the origin -- run a molecule for more");
    auto s=std::make_shared<AtomSetup>();
    s->opts=AtomOptions(asked);   // the kind's defaults filled, what an atom does not consume refused
    s->acc =acc;
    RunOptions& o=s->opts;
    s->Z =int(std::lround(atom.GetNuclearCharge()));
    const int charge=int(std::lround(atom.GetNetCharge()));
    s->Ne=s->Z-charge;
    s->structure=std::make_shared<Atom>(s->Z, double(charge), rvec3_t(0,0,0));   // at the origin, whatever the caller's position
    // THE RADIAL-ONLY ALIAS (V-XC-ONE-TERM step 5, controller 2026-11-14: deck backward compatibility).  Before AngularKind::None the atom
    // spelled its radial-only set Lebedev at degree 0, the 1-direction rule: the SAME points (one direction at 4 pi).  On the atom that
    // spelling still reads as it always did, now under its own name, so the admission rule and the reports see what the mesh IS.
    if (o.xcMesh->angular==qcMesh::AngularKind::Lebedev && o.xcMesh->angularDegree==0) o.xcMesh->angular=qcMesh::AngularKind::None;
    if (o.xcMesh->angular==qcMesh::AngularKind::None) o.xcMesh->angularDegree=0;   // not read by the rule: one value, one mesh ID
    s->xcMesh=*o.xcMesh;
    // The seed: an atom seeds from the core guess unless told otherwise (it never uses the molecular SAD superposition), and a stated
    // SAD reads the run's own atom: the GTH pseudo-atom at its Zion (the valence library), else all-electron.
    using ChargeDensity::SeedStrategy;
    if (o.seed==SeedStrategy::Default) o.seed=SeedStrategy::CoreGuess;
    if (o.atom.pseudopotential) s->seedKey=ChargeDensity::SeedKey("LDA", {{s->Z, s->PseudoValence()}});

    // THE ACCELERATOR: the options' type (its one home, V-FACADE-PARITY P2; Auto = DIIS, Null refused), unless the atom block's json escape
    // hatch NAMES one (its "type"), read through AcceleratorTypeByName: an unknown name is REFUSED -- until 2026-11-17 the atom ran DIIS on it
    // silently.  `directmin` is kept as an alias there: GDM stepping from the first iteration.  The two stated and different are refused.
    if (o.atom.accelerator.is_object() && o.atom.accelerator.contains("type"))
    {
        const std::string stated=o.atom.accelerator.value("type", std::string());
        s->directMin = stated=="directmin" || stated=="DirectMin";
        const auto t = s->directMin ? SCFAccelerators::Type::GDM : AcceleratorTypeByName(stated, "qchem::AtomCalculation");
        if (o.accelerator!=SCFAccelerators::Type::Auto && o.accelerator!=t)
            policy::Refuse("AcceleratorType", "accelerator "+std::string(NameOf(o.accelerator))+" and atom.accelerator type \""+stated+"\"",
                           "the two name different accelerators -- state one");
        s->accelerator=t;
    }
    else s->accelerator=MolecularAcceleratorType(o.accelerator, "qchem::AtomCalculation");
    // A GGA (stated model, or an atom.xc override): built, but the radial atom route does not sample the gradient yet (GGA stage 1).
    RefuseUnroutedGradient(*o.model, "an atom");
    if (o.atom.xc) RefuseUnroutedGradient(*o.atom.xc, "an atom (atom.xc)");
    const bool densityFock = o.atom.pseudopotential || o.atom.xc.has_value() || Hamiltonian::IsDFT(*o.model);
    if (densityFock) RefuseGradientConsumers(o.vxcFit, o.xcRho, NameOf(s->accelerator));   // GDM needs H = dE/dD (pin 32)
    // Atoms keep the proven Z-scaled DIIS engage gate unless the caller pins one (DFT too: \c DIISEngageGate's density-Fock rule is
    // the molecule's), over the molecular defaults.  The json escape hatch merges over these when the graph is built.
    const double eMax=DIISEngageGate(acc.diis.engageBelowFD.value_or(0.0), s->Z, /*the atom keeps the Z-scaled gate, DFT too*/false);
    s->accParams=ResolveAcceleratorParams(acc, MolecularAcceleratorDefaults(eMax));
    return s;
}

} // namespace qchem
