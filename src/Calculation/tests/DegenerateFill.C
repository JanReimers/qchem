// File: Calculation/tests/DegenerateFill.C  The facades' degenerate-level fill policy beside a direct minimiser (D-DEGEN-FILL follow-up,
// user 2026-11-08: "Many users will have no idea that GDM only works for integers. So we need to guide them."): a GDM run with the
// fill unstated takes the integer fill and the `choose` row teaches why; a GDM run that STATES Shared gets a `warn` row naming the
// fallback that will run instead.  One test per case and facade; the rows are filed before the first SCF iteration.
#include <gtest/gtest.h>
#include <functional>
#include <string>
#include <vector>
import qchem.DegenerateFillReport;
import qchem.MoleculeCalculation;
import qchem.SolidCalculation;
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.SCFParams;
import qchem.SCFAccelerator.Factory;
import qchem.SCFCap;
import qchem.Mesh;
import qchem.Reporting;
import qchem.Types;

using namespace qchem;
using Fill=SCFParams::DegenerateLevelFill;

namespace
{
//! Every `policy` row named \a name anywhere in \a j.
std::vector<report::json> PolicyRows(const report::json& j, const std::string& name)
{
    std::vector<report::json> out;
    std::function<void(const report::json&)> walk=[&](const report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.contains("name") && r["name"]==name) out.push_back(r);
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(j);
    return out;
}
const report::json& FiledRun(const std::string& prefix)
{
    const auto& G=report::GlobalReport();
    const report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(prefix,0)==0) last=&it.value();
    if (!last) throw std::logic_error("no filed run "+prefix);
    return *last;
}
bool Has(const report::json& row, const char* key, const std::string& text)
{
    return row.contains(key) && row[key].get<std::string>().find(text)!=std::string::npos;
}
SCFParams Bounded(Fill f)
{
    SCFParams p; p.NMaxIter=2; p.DegenerateLevels=f;
    p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
    return p;
}
//! Boron (2p^1) under LDA: the facade's own default is Shared -- what a stated GDM must override.
std::shared_ptr<const MoleculeSetup> BoronGDM()
{
    Molecule b; b.Insert(Site{5, Vector3D<double>(0,0,0)});
    return MoleculeSetup::Resolve(b, {.basis = "dzvp", .model = Model::LDA, .multiplicity = 2}, {.type = "GDM"});
}
//! Si Γ, small, stated grids (no pilot), GDM.
std::shared_ptr<const SolidSetup> SiGammaGDM(const char* label)
{
    const auto mat=Materials::Get("Si_diamond");
    SolidCalcOptions o; o.Nelec=8; o.densityEcut=20.0; o.label=label;
    o.accuracy.eps=0.0; o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    o.accelerator=SCFAccelerators::Type::GDM;
    return SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o);
}
} // namespace

TEST(DegenerateFillPolicy, TheDeclinedFallbackIsNamedWithTheRelaxItRunsAt)
{
    SCFParams p; p.StartingRelaxRo=0.5;
    p.Mixer=SCFParams::Mixing::Linear;
    EXPECT_EQ(DescribeDeclinedFallback(p), "linear D-mixing at relax 0.25") << "under the DeclinedDirectMinRelax ceiling";
    p.Mixer=SCFParams::Mixing::Pulay; p.KerkerG0=0.0; p.PulayDepth=8;
    EXPECT_EQ(DescribeDeclinedFallback(p), "Pulay D-mixing (depth 8) at relax 0.5") << "a history mixer is not capped (measured)";
    SCFParams a;    // the periodic Auto mixer: what a solid GDM run falls back to, at the relax the iterator resolves
    EXPECT_EQ(DescribeDeclinedFallback(ResolveMixing(a, /*gSpace*/true)), "Pulay D-mixing (depth 8) at relax 0.45");
}

TEST(DegenerateFillPolicy, AMoleculeWithGDMAndTheFillUnstatedTakesTheIntegerFillAndSaysWhy)
{
    {
        report::Run run("degen gdm molecule auto");
        BoundedSCF bounded("the policy rows are filed before the first iteration; the SCF is not the subject");
        MoleculeCalculation c(BoronGDM(), Bounded(Fill::Auto));
    }
    const auto rows=PolicyRows(FiledRun("degen gdm molecule auto@"), "DegenerateLevelFill");
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_TRUE(Has(rows[0], "decided", "IntegerOccupancy")) << rows[0].dump(1);
    EXPECT_TRUE(Has(rows[0], "reason", "GDM works only on integer occupations")) << rows[0].dump(1);
    EXPECT_TRUE(Has(rows[0], "reason", "degenerateLevelFill: Shared")) << "the row names how to get the symmetric state";
    EXPECT_TRUE(PolicyRows(FiledRun("degen gdm molecule auto@"), "DegenerateLevelFillUnderGDM").empty());
}

TEST(DegenerateFillPolicy, AMoleculeWithGDMAndSharedStatedIsWarnedWhatWillRun)
{
    {
        report::Run run("degen gdm molecule shared");
        BoundedSCF bounded("the policy rows are filed before the first iteration; the SCF is not the subject");
        MoleculeCalculation c(BoronGDM(), Bounded(Fill::Shared));
    }
    const auto& R=FiledRun("degen gdm molecule shared@");
    const auto choose=PolicyRows(R, "DegenerateLevelFill");
    ASSERT_EQ(choose.size(), 1u);
    EXPECT_TRUE(Has(choose[0], "decided", "SharedDegenerateOccupancy"));
    EXPECT_TRUE(Has(choose[0], "reason", "stated by the caller"));
    const auto warn=PolicyRows(R, "DegenerateLevelFillUnderGDM");
    ASSERT_EQ(warn.size(), 1u);
    EXPECT_EQ(warn[0]["verdict"], "warn");
    EXPECT_TRUE(Has(warn[0], "decided", "GDM will DECLINE")) << warn[0].dump(1);
    EXPECT_TRUE(Has(warn[0], "decided", "linear D-mixing at relax 0.25")) << "the molecule's fallback, at the relax it runs at";
}

TEST(DegenerateFillPolicy, ASolidWithGDMAndTheFillUnstatedTakesTheIntegerFillAndSaysWhy)
{
    {
        report::Run run("degen gdm solid auto");
        BoundedSCF bounded("the policy rows are filed before the first iteration; the SCF is not the subject");
        SolidCalculation c(SiGammaGDM("degen gdm solid auto"), Bounded(Fill::Auto));
    }
    const auto rows=PolicyRows(FiledRun("degen gdm solid auto@"), "DegenerateLevelFill");
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_TRUE(Has(rows[0], "decided", "IntegerOccupancy")) << rows[0].dump(1);
    EXPECT_TRUE(Has(rows[0], "reason", "GDM works only on integer occupations")) << rows[0].dump(1);
    EXPECT_TRUE(PolicyRows(FiledRun("degen gdm solid auto@"), "DegenerateLevelFillUnderGDM").empty());
}

TEST(DegenerateFillPolicy, ASolidWithGDMAndSharedStatedIsWarnedWhatWillRun)
{
    {
        report::Run run("degen gdm solid shared");
        BoundedSCF bounded("the policy rows are filed before the first iteration; the SCF is not the subject");
        SolidCalculation c(SiGammaGDM("degen gdm solid shared"), Bounded(Fill::Shared));
    }
    const auto& R=FiledRun("degen gdm solid shared@");
    const auto choose=PolicyRows(R, "DegenerateLevelFill");
    ASSERT_EQ(choose.size(), 1u);
    EXPECT_TRUE(Has(choose[0], "decided", "SharedDegenerateOccupancy"));
    const auto warn=PolicyRows(R, "DegenerateLevelFillUnderGDM");
    ASSERT_EQ(warn.size(), 1u);
    EXPECT_EQ(warn[0]["verdict"], "warn");
    EXPECT_TRUE(Has(warn[0], "decided", "GDM will DECLINE")) << warn[0].dump(1);
    EXPECT_TRUE(Has(warn[0], "decided", "instead")) << warn[0].dump(1);
}

//! Without GDM nothing changes: the facade's own default and no warning.
TEST(DegenerateFillPolicy, WithoutGDMTheFacadeDefaultStandsAndNothingWarns)
{
    report::Run run("degen nogdm pure");
    const SCFParams p=ResolveDegenerateLevelFill(SCFParams{}, Fill::Shared, "the facade's reason", true, {});
    EXPECT_EQ(p.DegenerateLevels, Fill::Shared);
    const SCFParams q=ResolveDegenerateLevelFill(SCFParams{}, Fill::Shared, "the facade's reason", true,
                                                 DegenerateFillContext{.directMin=true, .fallback="linear D-mixing at relax 0.25"});
    EXPECT_EQ(q.DegenerateLevels, Fill::OrbitalByOrbital) << "GDM overrides the facade's Shared default";
}
