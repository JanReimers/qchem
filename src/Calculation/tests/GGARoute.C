// File: Calculation/tests/GGARoute.C  GGA stage 1 (user ruling 2026-11-20): the GGA models are named, map to their libxc ids and to the
// GTH PPs fitted with them, and EVERY XC route refuses them by policy (GGARoute) until it samples the density gradient -- never an LDA in
// their place.
#include <gtest/gtest.h>
#include <filesystem>
#include <functional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
import qchem.MoleculeCalculation;
import qchem.AtomCalculation;
import qchem.SolidCalculation;
import qchem.Hamiltonian.Factory;   // Model, ModelToXC, LibxcIDs, ModelForGTHFitted, ConsumesDensityGradient
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;   // BasisSetData::SIPP_SR
import qchem.Structure;
import qchem.Types;
import qchem.Mesh;                    // UnitCellKind, BeckeXCParams: the solid route families
import qchem.Reporting;               // the filed run's policy row
import qchem.Deck;                    // a deck stating model PBE

using namespace qchem;

TEST(GGARoute, TheGGAModelsAreTheirLibxcPairsAndMatchTheirGTHFits)
{
    using H=Hamiltonian::Model;
    for (const auto& [m, name, ids] : std::vector<std::tuple<H,std::string,std::vector<int>>>{
             {H::PBE, "PBE", {101,130}}, {H::PBEsol, "PBESOL", {116,133}}, {H::BLYP, "BLYP", {106,131}}})
    {
        EXPECT_TRUE(Hamiltonian::IsDFT(m));
        EXPECT_EQ(Hamiltonian::ModelFromName(Hamiltonian::Name(m)), m) << "a deck spells it back";
        const auto xc=Hamiltonian::ModelToXC(m);
        EXPECT_EQ(Hamiltonian::LibxcIDs(xc), ids);
        EXPECT_EQ(Hamiltonian::GTHFittedName(xc), name);
        EXPECT_EQ(Hamiltonian::ModelForGTHFitted(name), m) << "a GTH-" << name << " PP now names the model that matches it";
        EXPECT_TRUE(Hamiltonian::ConsumesDensityGradient(xc));
    }
    for (H m : {H::LDA, H::PADE, H::Xalpha}) EXPECT_FALSE(Hamiltonian::ConsumesDensityGradient(Hamiltonian::ModelToXC(m)));
}

namespace
{
//! The filed run's GGARoute policy rows (the report's \c policy arrays, any depth).
std::vector<report::json> GGARouteRows(const std::string& runName)
{
    const auto& G=report::GlobalReport();
    const report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(runName+"@",0)==0) last=&it.value();
    if (!last) throw std::logic_error("no filed run "+runName);
    std::vector<report::json> out;
    std::function<void(const report::json&)> walk=[&](const report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.contains("name") && r["name"]=="GGARoute") out.push_back(r);
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(*last);
    return out;
}

//! One XC route family, and the setup call that would take a GGA down it.
struct RouteCase
{
    std::string label;
    std::function<std::string()> admit;   //!< the refusal's policy name ("" if it ADMITTED)
};

template <class O> std::string PolicyOf(const O& o) {return o ? std::string() : o.Error().policy;}

RunOptions SiPBE()
{
    RunOptions o; o.Nelec=8; o.solid.densityEcut=20.0; o.model=Hamiltonian::Model::PBE;
    return o;
}
std::string AdmitSi(const RunOptions& o)
{
    const auto si=Materials::Get("Si_diamond");
    return PolicyOf(SolidSetup::Admit(Lattice_3D(*si.cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o));
}
std::string AdmitHe(Hamiltonian::Model m, Hamiltonian::VxcFit fit)
{
    Molecule he; he.Insert(Site{2, rvec3_t(0,0,0)});
    RunOptions o; o.model=m; o.vxcFit=fit; o.molecule.basis="dzvp";
    return PolicyOf(MoleculeSetup::Admit(he, o));
}

std::vector<RouteCase> RouteFamilies()
{
    using H=Hamiltonian::Model; using F=Hamiltonian::VxcFit;
    std::vector<RouteCase> cs;
    for (H m : {H::PBE, H::PBEsol, H::BLYP})
    {
        const std::string n=Hamiltonian::Name(m);
        cs.push_back({"molecule Becke (" +n+")", [m]{return AdmitHe(m, F::Auto);}});
        cs.push_back({"molecule delta ("  +n+")", [m]{return AdmitHe(m, F::Delta);}});
        cs.push_back({"molecule auxiliary ("+n+")", [m]{return AdmitHe(m, F::Auxiliary);}});
        cs.push_back({"atom radial ("     +n+")", [m]{RunOptions o; o.model=m; return PolicyOf(AtomSetup::Admit(Atom(8), o));}});
        cs.push_back({"atom atom.xc ("    +n+")", [m]{RunOptions o; o.atom.xc=Hamiltonian::ModelToXC(m); return PolicyOf(AtomSetup::Admit(Atom(8), o));}});
    }
    cs.push_back({"GPW default mesh (Auto)", []{return AdmitSi(SiPBE());}});
    cs.push_back({"GPW raster (Uniform, plane-wave fit)", []{RunOptions o=SiPBE(); o.xcMesh=qcMesh::XCMeshDefault(); o.xcMesh->cellKind=qcMesh::UnitCellKind::Uniform; o.vxcFit=Hamiltonian::VxcFit::PlaneWave; return AdmitSi(o);}});
    cs.push_back({"GPW Becke solid (delta fit)", []{RunOptions o=SiPBE(); o.xcMesh=qcMesh::BeckeXCParams(); o.vxcFit=Hamiltonian::VxcFit::Delta; return AdmitSi(o);}});
    cs.push_back({"plane-wave orbitals", []{RunOptions o=SiPBE(); o.solid.orbitalBasis=RunOptions::Solid::OrbitalBasis::PlaneWave; o.solid.planeWaveEcut=4.0; return AdmitSi(o);}});
    cs.push_back({"plane-wave orbitals, polarized", []{RunOptions o=SiPBE(); o.solid.orbitalBasis=RunOptions::Solid::OrbitalBasis::PlaneWave; o.solid.planeWaveEcut=4.0; o.spin=SpinGroup::Polarized; return AdmitSi(o);}});
    return cs;
}
} // namespace

// THE NO-SILENT-SUBSTITUTION INVARIANT (GGA stage 1): a GGA is REFUSED by name -- policy GGARoute, a "refuse" row in the run report's
// policy scope -- on EVERY XC route family, since none has a density-gradient sampler yet (stages 2-4).  This is the factory side's ONE
// run-time check (user 2026-11-20: the gradient is a TYPE): past it, a GGA term cannot be constructed without a GGA functional and a
// gradient sampler, and a GGA functional has no rho-only face to be misused through (GGA_XC.C's static_asserts) -- so there is no
// "ran with v_sigma dropped" state left to guard at run time.
TEST(GGARoute, EveryRouteRefusesAGGAUntilItSamplesTheGradient)
{
    for (const RouteCase& c : RouteFamilies())
    {
        const std::string run="ggaroute "+c.label;
        std::string policy;
        {
            report::Run r(run);
            policy=c.admit();
        }
        EXPECT_EQ(policy, "GGARoute") << c.label << ": refused by name (empty = ADMITTED, a silent substitution)";
        const auto rows=GGARouteRows(run);
        ASSERT_EQ(rows.size(), 1u) << c.label << ": one GGARoute row in the report's policy scope";
        EXPECT_EQ(rows[0]["verdict"], "refuse") << c.label;
        const std::string decided=rows[0]["decided"], reason=rows[0]["reason"];
        EXPECT_NE(decided.find("libxc"), std::string::npos) << c.label << ": the row names the functional: " << decided;
        EXPECT_NE(reason.find("state an LDA"), std::string::npos) << c.label << ": and says what to do: " << reason;
    }
}

// A DECK stating model PBE is refused the same way: by deck::Vet (vet-decks, before anything is built) and by deck::Prepare (what Run
// builds through), both naming the policy.
TEST(GGARoute, ADeckStatingPBEIsRefusedByName)
{
    deck::RunSpec spec;
    deck::FromJson(report::json{{"structure","Si_diamond"},{"model","PBE"},{"basis",{{"data","SIPP_SR"}}}}, spec);
    for (const auto& [door, call] : std::vector<std::pair<std::string,std::function<void()>>>{
             {"Vet",     [&]{(void)deck::Vet(spec, std::filesystem::temp_directory_path());}},
             {"Prepare", [&]{(void)deck::Prepare(spec);}}})
    {
        try
        {
            call();
            ADD_FAILURE() << "deck::" << door << " ADMITTED a PBE deck: no route samples the gradient yet";
        }
        catch (const std::exception& e)
        {
            EXPECT_NE(std::string(e.what()).find("density gradient"), std::string::npos) << door << ": " << e.what();
        }
    }
}
