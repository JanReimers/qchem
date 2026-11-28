// File: Calculation/tests/AuditFootprint.C  The audit RSS policy (OpenWork (e), ruled 2026-11-07): the byte count, the calibration that
// must hold its own measurement record, the three verdicts, the memory query, and the two doors a refusal comes out of (the pilot's
// preflight, before any SCF; a postSCF audit).
#include <gtest/gtest.h>
#include <cmath>
#include <functional>
#include <string>
#include <vector>
import qchem.AuditFootprint;
import qchem.SolidCalculation;
import qchem.Materials;
import qchem.Lattice_3D;
import qchem.UnitCell;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.SCFParams;
import qchem.Mesh;
import qchem.Reporting;
import qchem.Types;

using namespace qchem;

TEST(AuditFootprint, TheByteCountIsTheTablesAndFieldsPerBlock)
{
    AuditReferenceSize s;
    s.rasterPoints=1000.0; s.meshPoints=200.0;
    s.blocks={{10,true},{12,false}};   // 8 B and 16 B an element
    const AuditBytes b=CountAuditBytes(s);
    EXPECT_DOUBLE_EQ(b.phiTables,    200.0*(10*8.0+12*16.0));
    EXPECT_DOUBLE_EQ(b.rasterFields, 1000.0*(8.0+16.0));
    s.spinChannels=2;
    EXPECT_DOUBLE_EQ(CountAuditBytes(s).rasterFields, 2.0*b.rasterFields) << "one raster field per block AND spin channel";
    EXPECT_DOUBLE_EQ(CountAuditBytes(s).phiTables, b.phiTables) << "the Phi tables are the basis's, not a spin channel's";
}

TEST(AuditFootprint, TheEstimateCombinesTheCountAndTheCalibrationInOnePlace)
{
    AuditReferenceSize s; s.rasterPoints=5.0e5; s.meshPoints=2.0e4; s.blocks={{23,true},{23,false}};
    const AuditFootprintCalibration cal{.eta=2.0, .kappa=100.0, .bandLow=0.5, .bandHigh=3.0};
    const AuditFootprint f=EstimateAuditFootprint(s, cal);
    EXPECT_DOUBLE_EQ(f.estimate, 2.0*CountAuditBytes(s).Known()+100.0*5.0e5);
    EXPECT_DOUBLE_EQ(f.Low(),  0.5*f.estimate);
    EXPECT_DOUBLE_EQ(f.High(), 3.0*f.estimate);
}

//! The calibration's own record (20 audits on 4 real runs) must sit inside its band: a refit that moves eta or kappa and leaves a measured
//! audit outside [bandLow, bandHigh] has made "refuse at the low end" unsafe, and fails here.
TEST(AuditFootprint, TheCalibrationHoldsItsOwnMeasurementRecord)
{
    const AuditFootprintCalibration& cal=kAuditFootprintCalibration;
    double lo=1e300, hi=0.0;
    for (const AuditFootprintMeasurement& m : kAuditFootprintRecord)
    {
        AuditReferenceSize s; s.rasterPoints=m.rasterPoints; s.meshPoints=m.meshPoints;
        for (int i=0;i<m.realBlocks;   ++i) s.blocks.push_back({size_t(m.nFunctions), true});
        for (int i=0;i<m.complexBlocks;++i) s.blocks.push_back({size_t(m.nFunctions), false});
        const double ratio=m.measuredMiB*1048576.0/EstimateAuditFootprint(s).estimate;
        EXPECT_GE(ratio, cal.bandLow)  << m.run << ", " << m.reference;
        EXPECT_LE(ratio, cal.bandHigh) << m.run << ", " << m.reference;
        lo=std::min(lo, ratio); hi=std::max(hi, ratio);
    }
    EXPECT_GT(lo, cal.bandLow-0.02)  << "the low band is the record's smallest ratio, rounded down -- not a wide margin";
    EXPECT_LT(hi, cal.bandHigh+0.02) << "the high band is the record's largest ratio, rounded up";
}

namespace
{
//! A room with a stated limit, own memory and live figure (bytes).
ResourceRoom Room(double limit, double own, double live)
{
    ResourceRoom r;
    r.limit={limit, "stated limit"};
    r.ownBytes=own; r.ownSource="stated own";
    r.live={live, "stated live"};
    return r;
}
}

TEST(AuditFootprint, TheThreeVerdicts)
{
    AuditFootprint f; f.estimate=100.0; f.cal={.eta=1.0, .kappa=0.0, .bandLow=0.8, .bandHigh=1.5};
    EXPECT_EQ(Judge(f, Room(200.0, 0.0, 200.0)), FootprintVerdict::Fits);
    EXPECT_EQ(Judge(f, Room(150.0, 0.0, 150.0)), FootprintVerdict::Fits)   << "the high end exactly fits";
    EXPECT_EQ(Judge(f, Room(120.0, 0.0, 120.0)), FootprintVerdict::Warn)   << "only the high end is over: it may not fit";
    EXPECT_EQ(Judge(f, Room( 80.0, 0.0,  80.0)), FootprintVerdict::Warn)   << "the low end exactly fits";
    EXPECT_EQ(Judge(f, Room( 79.0, 0.0,  79.0)), FootprintVerdict::Refuse) << "even the low end is over: it cannot fit";
}

// USER RULING 2026-11-27 (option B, "Fixed limit"): the refusal reads the deterministic limit less the run's own memory, never what the
// scope's neighbours and its page cache hold; the live figure only warns.  Lane H's panel (2026-11-26) refused NaF k222 relax, 0.40 GB,
// at "0.00 GB available" in a 9G scope.
TEST(AuditFootprint, TheRefusalReadsTheFixedRoomAndTheLiveFigureOnlyWarns)
{
    AuditFootprint f; f.estimate=100.0; f.cal={.eta=1.0, .kappa=0.0, .bandLow=0.8, .bandHigh=1.5};
    EXPECT_EQ(Judge(f, Room(200.0, 0.0, 0.0)), FootprintVerdict::Warn)    << "a full scope beside the run warns, it never refuses";
    EXPECT_EQ(Judge(f, Room(200.0, 0.0, 0.0)), Judge(f, Room(200.0, 0.0, 50.0))) << "the verdict does not move with the neighbours";
    EXPECT_EQ(Judge(f, Room(200.0, 121.0, 200.0)), FootprintVerdict::Refuse) << "the run's OWN memory counts: 79 left, the low end 80";
    EXPECT_EQ(Judge(f, Room(200.0, 40.0, 200.0)), FootprintVerdict::Fits)  << "160 left in the fixed room (the high end 150), the live figure clear";
    EXPECT_EQ(Judge(f, Room(200.0, 60.0, 200.0)), FootprintVerdict::Warn)  << "140 left: the high end 150 is over the fixed room";
    EXPECT_DOUBLE_EQ(Room(10.0, 20.0, 5.0).Fixed(), 0.0) << "the room is never negative";
    // Lane S, 2026-11-27: GPW_MnO.Γ_Shub_Pol_Smear_IonicSAD_eqSAD_Long under `memsafe -H 3G` alone in its scope, peak RSS 1.40 GB, was
    // refused its pilot's round-1 reference (0.66-1.27 GB) at "0.00 GB is available": the scope's page cache.  3 - 1.4 = 1.6 GB holds it.
    const double GB=1.0e9;
    AuditFootprint mno; mno.estimate=1.0*GB; mno.cal={.eta=1.0, .kappa=0.0, .bandLow=0.66, .bandHigh=1.27};
    EXPECT_EQ(Judge(mno, Room(3.0*GB, 1.40*GB, 0.0)), FootprintVerdict::Warn) << "a cap just above the run's own RSS does not refuse a run that fits";
    EXPECT_EQ(Judge(mno, Room(3.0*GB, 1.40*GB, 2.0*GB)), FootprintVerdict::Fits);
    const std::string d=DescribeRoom(Room(9e9, 1e9, 2e9));
    for (const char* s : {"fixed room 8.00 GB", "stated limit 9.00 GB", "own 1.00 GB", "live 2.00 GB", "stated live"})
        EXPECT_NE(d.find(s), std::string::npos) << s << " in: " << d;
}

TEST(AuditFootprint, TheMemoryFiguresParse)
{
    EXPECT_FALSE(ParseCgroupLimit("max\n"));
    EXPECT_DOUBLE_EQ(*ParseCgroupLimit("9663676416\n"), 9663676416.0);
    EXPECT_THROW((void)ParseCgroupLimit("lots\n"), std::runtime_error);
    EXPECT_THROW((void)ParseCgroupLimit(""), std::runtime_error);
    const std::string meminfo="MemTotal:       14000000 kB\nMemFree:          100000 kB\nMemAvailable:    6000000 kB\n";
    EXPECT_DOUBLE_EQ(*ParseMemAvailable(meminfo), 6000000.0*1024.0);
    EXPECT_FALSE(ParseMemAvailable("MemTotal: 1 kB\n"));
    // memory.stat counts bytes; only inactive_file and slab_reclaimable are reclaimable here (active_file stays used)
    const std::string stat="anon 1359953920\nfile 1466699776\nactive_file 241168384\ninactive_file 1224663040\nslab_reclaimable 101873400\n"
                           "slab_unreclaimable 5000000\n";
    EXPECT_DOUBLE_EQ(*ParseCgroupReclaimable(stat), 1224663040.0+101873400.0);
    EXPECT_FALSE(ParseCgroupReclaimable("anon 1\nfile 2\n"));
    const std::string status="Name:\tITMain\nVmRSS:\t  3000 kB\nRssAnon:\t  2000 kB\nRssFile:\t   900 kB\nRssShmem:\t   100 kB\n";
    EXPECT_DOUBLE_EQ(*ParseOwnResidentNonFile(status), 2100.0*1024.0) << "RssAnon + RssShmem; RssFile is page cache";
    EXPECT_FALSE(ParseOwnResidentNonFile("Name:\tx\n"));
}

TEST(AuditFootprint, TheMemoryAvailableNamesItsSource)
{
    const AvailableMemory a=QueryAvailableMemory();
    EXPECT_GE(a.bytes, 0.0);
    EXPECT_TRUE(a.source.rfind("cgroup memory.", 0)==0 || a.source=="MemAvailable (/proc/meminfo)") << a.source;
    const ResourceRoom r=QueryResourceRoom();
    EXPECT_GT(r.ownBytes, 0.0) << r.ownSource;
    EXPECT_GT(r.limit.bytes, r.ownBytes) << r.limit.source;
    EXPECT_GT(r.Fixed(), 0.0);
}

namespace
{
SCFParams Tight()
{
    SCFParams p;
    p.NMaxIter=100; p.MinΔρPerElectron=1e-8; p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
    return p;
}
//! Si Γ, small: a run whose pilot would run (accuracy > 0) or not (accuracy 0, the ball's own rasters stated), with a calibration so
//! large that ANY reference is refused -- the doors are gated, not the machine's memory.
std::shared_ptr<const SolidSetup> SiGamma(bool pilot, const SCFParams& scf)
{
    const auto mat=Materials::Get("Si_diamond");
    RunOptions o; o.Nelec=8; o.solid.densityEcut=20.0; o.label="auditMemory"; o.scf=scf;
    if (!pilot) { o.accuracy.eps=0.0; SolidXCMesh(o).cellKind=qcMesh::UnitCellKind::Uniform; }
    o.solid.auditFootprint.eta=1.0e9;   // ~1e9 x the byte count: no machine has it
    return SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=BasisSet::Gaussian::BasisSetData::SIPP_SR}, o);
}
//! Every `policy` row named \a name anywhere in a filed run's tree.
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
const report::json* FiledRun(const std::string& prefix)
{
    const auto& G=report::GlobalReport();
    const report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(prefix,0)==0) last=&it.value();
    return last;
}
} // namespace

//! A pilot that cannot fit is refused BEFORE the first SCF: no iteration is spent, the failure carries the numbers and the escape, and
//! the run report holds a `resource` policy row with the byte count, the calibration, the band and the room (both figures, their sources).
TEST(AuditFootprint, APilotThatCannotFitIsRefusedBeforeTheFirstSCF)
{
    {
        report::Run run("audit memory preflight");
        auto c=SolidCalculation::Prepare(SiGamma(/*pilot*/true, Tight())).TakeValue(); (void)c->Converge();
        const auto R=c->Result();
        ASSERT_FALSE(R);
        const SCFFailure& f=R.Error();
        EXPECT_EQ(f.why, SCFFailure::Why::ResourceRefused);
        EXPECT_EQ(f.iterations, 0u) << "refused before any iteration";
        ASSERT_TRUE(f.resource);
        EXPECT_EQ(f.resource->when, "before the first SCF");
        EXPECT_GT(f.resource->footprint.Low(), f.resource->room.Fixed()) << "refused on the fixed room";
        EXPECT_NE(f.resource->escape.find("accuracy 0"), std::string::npos) << "the refusal names the escape: " << f.resource->escape;
        EXPECT_FALSE(c->PilotCost().settled);
    }
    const auto* run=FiledRun("audit memory preflight@");
    ASSERT_NE(run, nullptr);
    const auto rows=PolicyRows(*run, "AuditReferenceMemory");
    ASSERT_EQ(rows.size(), 1u) << "the first reference refused; the second is never judged";
    const auto& r=rows.front();
    EXPECT_EQ(r["verdict"], "resource");
    for (const char* k : {"byteCount", "calibration", "band", "estimateBytes", "fixedRoomBytes", "limitSource", "ownBytes", "ownSource", "liveBytes",
                         "liveSource"})
        EXPECT_TRUE(r["detail"].contains(k)) << k << " missing from " << r.dump(1);
}

//! A postSCF audit that cannot fit is refused at its door (the run's own answer untouched); a run with stated grids has no pilot to refuse.
TEST(AuditFootprint, APostSCFAuditThatCannotFitIsRefusedAndTheRunStands)
{
    auto c=SolidCalculation::Prepare(SiGamma(/*pilot*/false, Tight())).TakeValue(); (void)c->Converge();
    ASSERT_TRUE(c->Result()) << "stated grids (accuracy 0): no pilot, nothing to refuse";
    const double E=c->Result()->Energy();
    const auto G=c->GridAudit();
    ASSERT_FALSE(G);
    EXPECT_EQ(G.Error().when, "postSCF gridAudit");
    const auto X=c->XCMeshAudit(qcMesh::MeshParams{.nUniform=24, .cellKind=qcMesh::UnitCellKind::Uniform});
    ASSERT_FALSE(X);
    EXPECT_EQ(X.Error().when, "postSCF xcMeshAudit");
    ASSERT_TRUE(c->Result());
    EXPECT_EQ(c->Result()->Energy(), E) << "a refused audit builds nothing and leaves the run's answer alone";
}
