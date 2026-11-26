// File: Calculation/tests/MoleculeSetup.C  The molecular RUN SETUP (doc/CalcNetworkReview.md §7, step 7): the decisions made
// once -- spin promotion, the valence electron count from the pseudopotential set, the parity and accelerator refusals.
#include <gtest/gtest.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <functional>
#include <nlohmann/json.hpp>
import qchem.MoleculeCalculation;
import qchem.Structure;
import qchem.Types;
import qchem.ChargeDensity.Seed;   // SeedStrategy, SeedLibrary (via AtomicDensity)
import qchem.SeedPolicy;           // SeedLibraryCovers: the Auto rule's coverage
import qchem.Vector3D;             // norm(rvec3_t)
import qchem.Structure.MeshFactory; // MakeMesh: what a point's maker must build (V-MESH-FACTORY)
import qchem.Reporting;             // report::Run / GlobalReport: the MeshSymmetry rows

using namespace qchem;

namespace
{
Molecule OneAtom(int Z)
{
    Molecule m;
    m.Insert(Site{Z, rvec3_t(0,0,0)});
    return m;
}
} // namespace

TEST(MoleculeSetup, AnOpenShellIsPromotedToPolarizedAndItsParityIsChecked)
{
    const auto s=MoleculeSetup::Resolve(OneAtom(7), {.multiplicity=4, .molecule={.basis="dzvp"}});
    EXPECT_EQ(s->opts.spin, SpinGroup::Polarized);
    EXPECT_EQ(s->Ne, 7);
    EXPECT_FALSE(s->pp) << "no pp named: all-electron";
    const auto refused=MoleculeSetup::Admit(OneAtom(7), {.multiplicity=3, .molecule={.basis="dzvp"}});
    ASSERT_FALSE(refused) << "a triplet of 7 electrons";
    EXPECT_EQ(refused.Error().policy, "SpinMultiplicity");
    EXPECT_THROW((void)MoleculeSetup::Resolve(OneAtom(7), {.multiplicity=3, .molecule={.basis="dzvp"}}), std::invalid_argument) << "the code door throws it";
}

TEST(MoleculeSetup, APseudopotentialSetGivesTheValenceCountAndTheValenceStructure)
{
    const Molecule na=OneAtom(11);
    const auto s=MoleculeSetup::Resolve(na, {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="valence_lowq_sr"}});
    ASSERT_TRUE(s->pp);
    EXPECT_EQ(s->Ne, s->pp->Valence("Na"));
    const auto point=s->MakePoint().structure;
    EXPECT_EQ(int(point->GetNumElectrons()), s->Ne);
    EXPECT_EQ(int(s->structure->GetNumElectrons()), 11) << "the reference keeps the TRUE charges; only a point is valence-charged";
}

// A CHARGED molecule keeps its charge under a pseudopotential (lane G, 2026-11-21): the point structure REPLACED each site's stated
// net charge with Z - Zion, so OH- with GTH ran NEUTRAL, 7 valence electrons instead of 8 -- a silent wrong answer.  A point now ADDS
// the core charge Z - Zion to the stated one: the valence count is sum(Zion) minus the stated charge, O q6 + H q1 + 1 = 8.
TEST(MoleculeSetup, AChargedMoleculeKeepsItsChargeUnderAPseudopotential)
{
    Molecule oh;
    oh.Insert(Site{8, rvec3_t(0,0,0), -1.0});        // OH-: the anion's extra electron on the O site
    oh.Insert(Site{1, rvec3_t(0,0,1.8330)});
    const auto s=MoleculeSetup::Resolve(oh, {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="dzvp_molopt_gth"}});
    ASSERT_TRUE(s->pp);
    const int valence=s->pp->Valence("O")+s->pp->Valence("H");
    EXPECT_EQ(s->Ne, valence+1) << "OH- carries one electron past its valence count";
    const auto point=s->MakePoint().structure;
    EXPECT_EQ(int(point->GetNumElectrons()), s->Ne);
    EXPECT_DOUBLE_EQ(point->GetNetCharge()-s->structure->GetNetCharge(), double(9-valence))
        << "a point adds the core charges Z - Zion to the stated ones";
    EXPECT_EQ(int(s->structure->GetNumElectrons()), 10) << "the reference keeps the TRUE charges: 9 + 1";
}

// A pseudopotential run seeds its OWN pseudo-atoms (OpenWork §3 "A pseudopotential run seeds from the VALENCE library for its q, on every
// structure", lane B 2026-11-12): the setup keys the seed library by each element's PP q, so a PP molecule reads the VALENCE densities at
// that q, as a solid does, and an all-electron one stays unkeyed (the all-electron library).  The tell: Si has no ALL-ELECTRON LDA entry,
// so a Si PP molecule used to take the core guess with a WARN; its q4 VALENCE entry exists, so it now takes SAD.
TEST(MoleculeSetup, APseudopotentialRunSeedsFromTheValenceLibraryAtItsQ)
{
    const Molecule si=OneAtom(14);
    const auto pp=MoleculeSetup::Resolve(si, {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="sipp"}});
    EXPECT_EQ(pp->seedKey.Q(14), pp->pp->Valence("Si")) << "the seed key carries the run's own q";
    EXPECT_EQ(ChargeDensity::SeedLibrary(pp->seedKey, 14), "atomic_valence_densities.json");
    EXPECT_EQ(pp->opts.seed, ChargeDensity::SeedStrategy::SAD) << "the valence library holds Si q4, so the Auto rule seeds from SAD";

    const auto ae=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::LDA, .molecule={.basis="dzvp"}});
    EXPECT_EQ(ae->seedKey.Q(8), -1) << "an all-electron run is unkeyed";
    EXPECT_EQ(ChargeDensity::SeedLibrary(ae->seedKey, 8), "atomic_densities.json");
}

// An ALL-ELECTRON molecule past neon seeds from SAD (OpenWork §3 "The default SAD seed stops at Ne", lane B 2026-11-13): the all-electron
// library held Z 1-10 only, so TeF6 on dzvp took the core guess with a WARN (before 3afed2228 it died on an uncaught throw).  It now holds
// every element NIST tabulates, so the Auto rule's coverage holds through uranium and the Te atom's setup picks SAD.
TEST(MoleculeSetup, AnAllElectronMoleculePastNeonSeedsFromSAD)
{
    const ChargeDensity::SeedKey unkeyed("LDA");
    for (int Z=1; Z<=92; Z++) EXPECT_TRUE(SeedLibraryCovers(unkeyed, Z)) << "Z=" << Z;
    const auto te=MoleculeSetup::Resolve(OneAtom(52), {.model=Model::LDA, .molecule={.basis="dzvp"}});
    EXPECT_EQ(ChargeDensity::SeedLibrary(te->seedKey, 52), "atomic_densities.json") << "an all-electron run reads the all-electron library";
    EXPECT_EQ(te->opts.seed, ChargeDensity::SeedStrategy::SAD) << "the library holds Te, so the Auto rule seeds from SAD";
}

TEST(MoleculeSetup, ABasisWithoutAnElementIsRefusedAtResolveNamingTheBasesThatCoverIt)
{
    try
    {
        (void)MoleculeSetup::Resolve(OneAtom(6), {.molecule={.basis="sipp"}});
        FAIL() << "sipp has no carbon: the setup must refuse it before anything is built";
    }
    catch (const std::invalid_argument& e)
    {
        const std::string what=e.what();
        EXPECT_NE(what.find("basis sipp has no block for C"), std::string::npos) << what;
        EXPECT_NE(what.find("dzvp"), std::string::npos) << "the registered bases that cover carbon are named: " << what;
    }
}

// The accelerator is a TYPE in its one home (RunOptions::accelerator, V-FACADE-PARITY P2): an unknown NAME is refused where text becomes a
// type (AcceleratorTypeByName: a GUI form, the atom's json escape hatch; the deck's enum reader), and the type a molecule cannot run, Null,
// at Admit.
TEST(MoleculeSetup, AnAcceleratorAMoleculeCannotRunIsRefusedAtResolve)
{
    const auto refused=MoleculeSetup::Admit(OneAtom(2), {.accelerator=SCFAccelerators::Type::Null, .molecule={.basis="dzvp"}});
    ASSERT_FALSE(refused);
    EXPECT_EQ(refused.Error().policy, "AcceleratorType");
    EXPECT_NE(refused.Error().reason.find("Null"), std::string::npos) << refused.Error().reason;
    try { (void)AcceleratorTypeByName("Anderson", "the test"); ADD_FAILURE() << "an unknown name became a type"; }
    catch (const std::exception& e) { EXPECT_NE(std::string(e.what()).find("Anderson"), std::string::npos) << e.what(); }
}

// EVERY KNOB OF THE MOLECULAR CHAIN REACHES THE FACTORY'S DOOR (2026-11-02; the typed door since V-FACADE-PARITY P0b-1): a molecule that
// states nothing runs what it always ran -- the molecular defaults over the Z-scaled gate -- and every stated knob lands in its field.
TEST(MoleculeSetup, EveryAcceleratorKnobReachesTheFactoryDoor)
{
    const SCFAccelerators::AcceleratorParams def=ResolveAcceleratorParams(AcceleratorOptions{}, MolecularAcceleratorDefaults(0.3));
    EXPECT_EQ(def.diis.Nproj, 4u); EXPECT_EQ(def.diis.EngageBelowFD, 0.3); EXPECT_EQ(def.diis.SVTol, 0.0); EXPECT_EQ(def.diis.SVTolRel, 1e-4);
    EXPECT_EQ(def.gdm.EngageBelowFD, 0.3) << "the standalone GDM gate follows the DIIS gate, as the flat eMax always did";
    EXPECT_EQ(def.gdm.Trust, 0.1); EXPECT_EQ(def.ladderGdm.EngageBelowFD, 1e10) << "the Ladder's GDM rung always steps";
    EXPECT_EQ(def.ladderEThresh, 1e-8); EXPECT_EQ(def.ladderStall, 5); EXPECT_EQ(def.ladderFloor, 1e-8); EXPECT_EQ(def.ladderSwitchAt, 1e-6);
    EXPECT_EQ(def.ladderSignal, SCFAccelerators::ScheduleSignal::EnergyChange);
    AcceleratorOptions a; a.diis.nProj=6; a.ladder.gdm.trust=0.25; a.ladder.eThresh=2e-8; a.ladder.stall=7; a.ladder.floor=3e-8; a.ladder.switchAt=1e-6;
    const SCFAccelerators::AcceleratorParams p=ResolveAcceleratorParams(a, MolecularAcceleratorDefaults(0.5));
    EXPECT_EQ(p.diis.Nproj, 6u); EXPECT_EQ(p.diis.EngageBelowFD, 0.5); EXPECT_EQ(p.ladderGdm.Trust, 0.25); EXPECT_EQ(p.gdm.Trust, 0.1) << "stated for the rung only";
    EXPECT_EQ(p.ladderEThresh, 2e-8); EXPECT_EQ(p.ladderStall, 7); EXPECT_EQ(p.ladderFloor, 3e-8); EXPECT_EQ(p.ladderSwitchAt, 1e-6);
}

// "EMax" was the [F,D] engage gate, never an energy: RETIRED from the accelerator json door (user ruling 2026-11-17), refused BY NAME with
// the new spelling in the message; the door's keys (ToAcceleratorJson) spell it EngageBelowFD.
TEST(MoleculeSetup, TheRetiredEMaxKeyIsRefusedNamingEngageBelowFD)
{
    using SCFAccelerators::Type;
    try { delete SCFAccelerators::Factory(Type::DIIS, nlohmann::json{{"NProj",4},{"EMax",0.1},{"SVTol",0.0}}); FAIL() << "EMax must be refused"; }
    catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find("EngageBelowFD"), std::string::npos) << e.what(); }
    const nlohmann::json j=ToAcceleratorJson(MolecularAcceleratorDefaults(0.2));
    EXPECT_FALSE(j.contains("EMax")); EXPECT_EQ(j.at("EngageBelowFD"), 0.2);
    std::unique_ptr<SCFAccelerators::SCFAccelerator> ok(SCFAccelerators::Factory(Type::DIIS, j));
    EXPECT_TRUE(ok);
}

// THE COULOMB FIT IS A SETUP POLICY (policy CoulombFit; OpenWork "The Coulomb-fit default is a hidden policy", lane S 2026-11-21): a density
// model's unstated Hartree route resolves HERE to the Auxiliary (Dunlap) fit on a NAMED file, and every use the name cannot have is refused
// by name -- before the row, the orbital families read A1_coul.bsd themselves and nothing above them could see or change it.
TEST(MoleculeSetup, AnUnstatedCoulombFitIsNamedByTheSetup)
{
    using CF=Hamiltonian::CoulombFit;
    const auto s=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::LDA, .molecule={.basis="dzvp"}});
    EXPECT_EQ(s->opts.coulomb.fit, CF::Delta) << "the molecular default (user ruling 2026-11-22): the route that meets eps";
    EXPECT_TRUE(s->opts.coulomb.auxBasis.empty());
    const auto aux=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::LDA, .coulomb={.fit=CF::Auxiliary}, .molecule={.basis="dzvp"}});
    EXPECT_EQ(aux->opts.coulomb.auxBasis, "A2_coul") << "Auxiliary stated without a file: A2, the best Gaussian fit measured";
    const auto te=MoleculeSetup::Resolve(OneAtom(52), {.model=Model::LDA, .coulomb={.fit=CF::Auxiliary}, .molecule={.basis="dzvp"}});
    EXPECT_EQ(te->opts.coulomb.auxBasis, "A1_coul") << "A2 stops at Zn: the unnamed Auxiliary fit falls back to A1 (H-Te), said so";
    const auto a2=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::LDA, .coulomb={.auxBasis="A2_coul"}, .molecule={.basis="dzvp"}});
    EXPECT_EQ(a2->opts.coulomb.fit, CF::Auxiliary) << "a stated fit file implies its route";
    EXPECT_EQ(a2->opts.coulomb.auxBasis, "A2_coul");
    const auto hf=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::HF, .molecule={.basis="dzvp"}});
    EXPECT_TRUE(hf->opts.coulomb.auxBasis.empty()) << "HF reads no fit: nothing is named";
    const auto delta=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::LDA, .coulomb={.fit=CF::Delta}, .molecule={.basis="dzvp"}});
    EXPECT_EQ(delta->opts.coulomb.fit, CF::Delta);
    EXPECT_TRUE(delta->opts.coulomb.auxBasis.empty());
    // V-FIT-AXES: an unstated axis on the Auxiliary route is Dunlap's (Coulomb metric, charge-constrained); a stated axis implies that route
    EXPECT_EQ(aux->opts.coulomb.metric, Hamiltonian::CoulombMetric::Coulomb);
    EXPECT_EQ(aux->opts.coulomb.constraint, Hamiltonian::ChargeConstraint::Charge);
    const auto plain=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::LDA, .coulomb={.constraint=Hamiltonian::ChargeConstraint::None}, .molecule={.basis="dzvp"}});
    EXPECT_EQ(plain->opts.coulomb.fit, CF::Auxiliary) << "a stated fit axis implies the Gaussian route, as a stated file does";
    EXPECT_EQ(plain->opts.coulomb.auxBasis, "A2_coul");
    EXPECT_EQ(plain->opts.coulomb.metric, Hamiltonian::CoulombMetric::Coulomb);
    EXPECT_EQ(plain->opts.coulomb.constraint, Hamiltonian::ChargeConstraint::None);
    EXPECT_EQ(delta->opts.coulomb.constraint, Hamiltonian::ChargeConstraint::Auto) << "the δ route has no constraint axis: left unresolved";
}

TEST(MoleculeSetup, ACoulombFitFileTheRunCannotUseIsRefusedByName)
{
    using CF=Hamiltonian::CoulombFit;
    struct Case { RunOptions o; std::string policy; int Z; };
    const std::vector<Case> cases={
        {{.model=Model::LDA, .coulomb={.fit=CF::Delta, .auxBasis="A1_coul"}, .molecule={.basis="dzvp"}}, "CoulombAuxBasisUnused", 8},
        {{.model=Model::LDA, .coulomb={.fit=CF::Multipole, .auxBasis="A1_coul"}, .molecule={.basis="dzvp"}}, "CoulombAuxBasisUnused", 8},
        {{.model=Model::HF, .coulomb={.fit=CF::Delta}, .molecule={.basis="dzvp"}}, "CoulombRouteForHF", 8},      // RI-J HF is the Gaussian fit's
        {{.model=Model::HF, .coulomb={.fit=CF::Multipole}, .molecule={.basis="dzvp"}}, "CoulombRouteForHF", 8},
        {{.model=Model::LDA, .coulomb={.auxBasis="no_such_fit"}, .molecule={.basis="dzvp"}}, "CoulombAuxBasisUnknown", 8},
        {{.model=Model::LDA, .coulomb={.auxBasis="A2_coul"}, .molecule={.basis="dzvp"}}, "CoulombAuxBasisCoverage", 52},   // A2 stops at Zn
        // V-FIT-AXES: the δ fit refuses the axes (user ruling 2026-11-25), the overlap metric is not built yet, HF fits nothing
        {{.model=Model::LDA, .coulomb={.fit=CF::Delta, .constraint=Hamiltonian::ChargeConstraint::Charge}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        {{.model=Model::LDA, .coulomb={.fit=CF::Multipole, .constraint=Hamiltonian::ChargeConstraint::Charge}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        {{.model=Model::LDA, .coulomb={.fit=CF::Delta, .metric=Hamiltonian::CoulombMetric::Coulomb}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        {{.model=Model::LDA, .coulomb={.fit=CF::Auxiliary, .metric=Hamiltonian::CoulombMetric::Overlap}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        // The analytic (ij|δ) solver: the δ fit's only (lane S 2026-11-25)
        {{.model=Model::HF, .coulomb={.solver=Hamiltonian::CoulombSolver::Analytic}, .molecule={.basis="dzvp"}}, "CoulombSolverUnused", 8},
        {{.model=Model::LDA, .coulomb={.fit=CF::Auxiliary, .solver=Hamiltonian::CoulombSolver::Analytic}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        {{.model=Model::LDA, .coulomb={.fit=CF::Multipole, .solver=Hamiltonian::CoulombSolver::Analytic}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        // The exact four-centre J (coulomb.route Exact): every fit key beside it is refused by name; HF states nothing (its J is exact)
        {{.model=Model::LDA, .coulomb={.route=Hamiltonian::HartreeRoute::Exact, .fit=CF::Delta}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        {{.model=Model::LDA, .coulomb={.route=Hamiltonian::HartreeRoute::Exact, .auxBasis="A2_coul"}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        {{.model=Model::LDA, .coulomb={.route=Hamiltonian::HartreeRoute::Exact, .constraint=Hamiltonian::ChargeConstraint::None}, .molecule={.basis="dzvp"}}, "AdmitCoulomb", 8},
        {{.model=Model::HF, .coulomb={.route=Hamiltonian::HartreeRoute::Exact}, .molecule={.basis="dzvp"}}, "CoulombRouteUnused", 8},
    };
    for (const Case& c : cases)
    {
        const auto r=MoleculeSetup::Admit(OneAtom(c.Z), c.o);
        ASSERT_FALSE(r) << c.policy;
        EXPECT_EQ(r.Error().policy, c.policy);
    }
}

// RI-J HF (V-HAM-COMPOSE step 2): HF's J is the exact four-centre one unless a Gaussian fit is stated (the fit, a file or a fit axis); then
// it resolves as a density model's Auxiliary fit, and its row names it.
TEST(MoleculeSetup, HartreeFockTakesTheExactJUnlessAGaussianFitIsStated)
{
    using CF=Hamiltonian::CoulombFit;
    const auto exact=MoleculeSetup::Resolve(OneAtom(8), {.model=Model::HF, .molecule={.basis="dzvp"}});
    EXPECT_EQ(exact->opts.coulomb.route, Hamiltonian::HartreeRoute::Exact);
    EXPECT_FALSE(Hamiltonian::UsesCoulombFitBasis(exact->opts.coulomb)) << "the exact J reads no fit basis";
    for (const RunOptions& o : {RunOptions{.model=Model::HF, .coulomb={.fit=CF::Auxiliary}, .molecule={.basis="dzvp"}},
                                RunOptions{.model=Model::HF, .coulomb={.auxBasis="A1_coul"}, .molecule={.basis="dzvp"}},
                                RunOptions{.model=Model::HF, .coulomb={.constraint=Hamiltonian::ChargeConstraint::None}, .molecule={.basis="dzvp"}}})
    {
        const auto rij=MoleculeSetup::Resolve(OneAtom(8), o);
        EXPECT_EQ(rij->opts.coulomb.route, Hamiltonian::HartreeRoute::Fit);
        EXPECT_EQ(rij->opts.coulomb.fit, CF::Auxiliary);
        EXPECT_FALSE(rij->opts.coulomb.auxBasis.empty()) << "a fit file resolved";
        EXPECT_EQ(rij->opts.coulomb.metric, Hamiltonian::CoulombMetric::Coulomb);
        EXPECT_TRUE(Hamiltonian::UsesCoulombFitBasis(rij->opts.coulomb));
    }
}

// A GEOMETRY POINT (V-MESH-FACTORY §8 step 2): MakePoint(R) moves the reference's sites to R, keeps its charges (valence-charged under a
// pseudopotential), and its mesh maker builds over THAT point -- never over the reference -- exactly the factory's mesh for it.
TEST(MoleculeSetup, APointsMeshMakerBuildsOverThePointNotTheReference)
{
    Molecule h2;
    h2.Insert(Site{1, rvec3_t(0,0,0)});
    h2.Insert(Site{1, rvec3_t(0,0,1.4)});
    const auto s=MoleculeSetup::Resolve(h2, {.model=Model::LDA, .molecule={.basis="dzvp"}});
    const std::vector<rvec3_t> R={rvec3_t(0.1,0,0), rvec3_t(0,0.2,1.6)};
    const MoleculeSetup::Point p=s->MakePoint(R);
    ASSERT_TRUE(p.structure); ASSERT_TRUE(p.meshMaker);
    auto positions=[](const Structure& st){ std::vector<rvec3_t> r; st.ForEachSite([&](int, const rvec3_t& x, bool){ r.push_back(x); }); return r; };
    const std::vector<rvec3_t> at=positions(*p.structure), ref=positions(*s->structure);
    for (size_t a=0; a<R.size(); a++) EXPECT_EQ(norm(at[a]-R[a]), 0.0) << "the point is at R";
    EXPECT_EQ(norm(ref[1]-rvec3_t(0,0,1.4)), 0.0) << "the reference does not move";
    qcMesh::MeshParams mp=s->xcMesh; mp.nRadial=20; mp.angularDegree=7;
    const qcMesh::Mesh got=p.meshMaker->Make(mp), want=MakeMesh(Molecule(*p.structure), mp);
    ASSERT_EQ(got.size(), want.size());
    ASSERT_EQ(got.NSites(), size_t(2));
    for (size_t q=0; q<got.size(); q++) ASSERT_EQ(norm(got.Points()[q]-want.Points()[q]), 0.0);
    const qcMesh::Mesh atRef=MakeMesh(h2, mp);
    EXPECT_GT(norm(got.Points()[0]-atRef.Points()[0]), 0.05) << "not the reference's mesh";
}

// POLICY MeshSymmetry (OOD-SOLID V-FOLD-NEUTRAL; user ruling 2026-11-23, D1-D3): a molecule's XC mesh is made invariant under the
// D2h-type subgroup of its point group (rotated into the group's frame), folded only when the run imposes symmetry; a C1 molecule, an
// angular rule that cannot carry the frame, and a geometry symmetric only to the detection tolerance keep the lab frame, each REPORTED.
// The source is per geometry POINT: a step that breaks a mirror gets the smaller group there, and says so.
namespace
{
std::vector<report::json> PolicyRows(const std::string& run, const std::string& name)
{
    std::vector<report::json> out;
    const auto& G=report::GlobalReport();
    const report::json* filed=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(run+"@",0)==0) filed=&it.value();
    if (!filed) throw std::logic_error(run+" was not filed");
    std::function<void(const report::json&)> walk=[&](const report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.value("name","")==name) out.push_back(r);
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(*filed);
    return out;
}
Molecule Water(double dy=0.0)
{
    Molecule h2o;
    h2o.Insert(Site{8, rvec3_t(0,0,0)});
    h2o.Insert(Site{1, rvec3_t(0, 1.43+dy,1.11)});
    h2o.Insert(Site{1, rvec3_t(0,-1.43,1.11)});
    return h2o;
}
} //anon

TEST(MoleculeSetup, AMoleculesMeshIsInvariantUnderItsD2hTypeSubgroupAndSaysSo)
{
    std::shared_ptr<const MoleculeSetup> lda, imposed, hf;
    { report::Run run("MeshSymmetry LDA");     lda    =MoleculeSetup::Resolve(Water(), {.model=Model::LDA, .molecule={.basis="dzvp"}}); }
    { report::Run run("MeshSymmetry imposed"); imposed=MoleculeSetup::Resolve(Water(), {.model=Model::LDA, .imposeSymmetry=true, .molecule={.basis="dzvp"}}); }
    { report::Run run("MeshSymmetry HF");      hf     =MoleculeSetup::Resolve(Water(), {.model=Model::HF,  .molecule={.basis="dzvp"}}); }
    const MoleculeSetup::Point p=lda->MakePoint(), q=imposed->MakePoint();
    ASSERT_TRUE(p.meshSymmetry && q.meshSymmetry && hf->MakePoint().meshSymmetry);
    EXPECT_EQ(p.meshSymmetry->MeshOps().size(), size_t(4)) << "C2v";
    EXPECT_FALSE(p.meshSymmetry->AveragesDensityOverOrbits()) << "a free run: invariant, not folded (D2)";
    EXPECT_EQ(q.meshSymmetry->MeshOps().size(), size_t(4)) << "an imposed run: the same invariant mesh";
    for (const auto& op : p.meshSymmetry->MeshOps()) EXPECT_EQ(norm(op.tau), 0.0) << "point-group ops act about the centroid the maker supplies";
    const auto r=PolicyRows("MeshSymmetry LDA", "MeshSymmetry");
    ASSERT_EQ(r.size(), size_t(1));
    EXPECT_NE(r[0].dump().find("C2v"), std::string::npos) << r[0].dump();
    EXPECT_NE(r[0].dump().find("not folded"), std::string::npos) << r[0].dump();
    EXPECT_NE(PolicyRows("MeshSymmetry imposed", "MeshSymmetry")[0].dump().find("FOLDED"), std::string::npos)
        << "an imposed δ-fit DFT run samples XC once per orbit (OpenWork §3, the per-orbit consumer)";
    EXPECT_TRUE(imposed->MakePoint().meshSymmetry->AveragesDensityOverOrbits());
    std::shared_ptr<const MoleculeSetup> a1;
    { report::Run run("MeshSymmetry imposed A1"); a1=MoleculeSetup::Resolve(Water(), {.model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary,
                                                                                       .imposeSymmetry=true, .molecule={.basis="dzvp"}}); }
    EXPECT_FALSE(a1->MakePoint().meshSymmetry->AveragesDensityOverOrbits()) << "a Gaussian Vxc fit's functions are not invariant: not folded";
    EXPECT_NE(PolicyRows("MeshSymmetry imposed A1", "MeshSymmetry")[0].dump().find("NOT folded"), std::string::npos);
    EXPECT_TRUE(PolicyRows("MeshSymmetry HF", "MeshSymmetry").empty()) << "an all-electron HF run builds no XC mesh";
    EXPECT_TRUE(hf->MakePoint().meshSymmetry->MeshOps().empty());
}

TEST(MoleculeSetup, AMoleculesMeshKeepsTheLabFrameWhereItMustAndSaysWhy)
{
    std::shared_ptr<const MoleculeSetup> rough, coarse;
    { report::Run run("MeshSymmetry rough"); rough=MoleculeSetup::Resolve(Water(1e-6), {.model=Model::LDA, .molecule={.basis="dzvp"}}); }
    EXPECT_TRUE(rough->MakePoint().meshSymmetry->MeshOps().empty()) << "symmetric to the detection tolerance only";
    const auto r=PolicyRows("MeshSymmetry rough", "MeshSymmetry");
    ASSERT_EQ(r.size(), size_t(1));
    EXPECT_NE(r[0].dump().find("lab frame"), std::string::npos) << r[0].dump();
    EXPECT_NE(r[0].dump().find("symmetric to more digits"), std::string::npos) << r[0].dump();

    RunOptions o{.model=Model::LDA, .molecule={.basis="dzvp"}};
    o.xcMesh=qcMesh::MeshParams{}; o.xcMesh->nRadial=40; o.xcMesh->angular=qcMesh::AngularKind::Lebedev; o.xcMesh->angularDegree=7;
    { report::Run run("MeshSymmetry coarse"); coarse=MoleculeSetup::Resolve(Water(), o); }
    EXPECT_TRUE(coarse->MakePoint().meshSymmetry->MeshOps().empty()) << "the 24-direction rule cannot carry the frame";
    EXPECT_NE(PolicyRows("MeshSymmetry coarse", "MeshSymmetry")[0].dump().find("degree 9"), std::string::npos);

    // A ROTATED rule (BeckeXCParams() carries the cubic-site calibration's angRot) is a stated orientation: kept, in the lab frame.
    std::shared_ptr<const MoleculeSetup> rotated;
    o.xcMesh=qcMesh::BeckeXCParams(); o.xcMesh->angularDegree=41;
    ASSERT_NE(o.xcMesh->angRot, 0.0) << "this case needs the calibrated rotation";
    { report::Run run("MeshSymmetry rotated"); rotated=MoleculeSetup::Resolve(Water(), o); }
    EXPECT_TRUE(rotated->MakePoint().meshSymmetry->MeshOps().empty()) << "a rotated rule cannot carry the frame";
    EXPECT_NE(PolicyRows("MeshSymmetry rotated", "MeshSymmetry")[0].dump().find("angRot"), std::string::npos);

    // A geometry step that breaks C2v's mirrors (one H moved off) keeps only the molecular plane: Cs, reported as a change of group.
    std::shared_ptr<const MoleculeSetup> lda;
    { report::Run run("MeshSymmetry step"); lda=MoleculeSetup::Resolve(Water(), {.model=Model::LDA, .molecule={.basis="dzvp"}});
      const MoleculeSetup::Point moved=lda->MakePoint({rvec3_t(0,0,0), rvec3_t(0,1.60,1.11), rvec3_t(0,-1.43,1.11)});
      EXPECT_EQ(moved.meshSymmetry->MeshOps().size(), size_t(2)) << "Cs"; }
    const auto s=PolicyRows("MeshSymmetry step", "MeshSymmetry");
    ASSERT_EQ(s.size(), size_t(2)) << "the reference's row, and the step's change of group";
    EXPECT_NE(s[1].dump().find("this geometry point"), std::string::npos) << s[1].dump();
}
