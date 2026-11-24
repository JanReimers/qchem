// File: ChargeDensity/tests/BlockInstanceKeys.C  A fitter's per-block views are keyed by the BLOCK, not by its symmetry (lane H 2026-11-15).
//
// The basis-on-grid survey (controller ruling 2026-11-15, P1): every fitter, sampler and projector that holds a view per orbital block
// keyed it by the block's SPATIAL symmetry, so two blocks of equal symmetry (two basis sets, or one rebuilt at a new geometry, against a
// surviving fitter) got the FIRST block's view -- the latent defect S3c removed from the δ fit basis's own Φ tables
// (M_BasisCache.TwoBlocksOfEqualSymmetryAskedOfOneDeltaFitGetTheirOwnTables).  No production path reaches it (one fitter per Calculation,
// and every geometry point or raster rebuild builds a new one), so this gate is at the library level, in S3c's shape: two waters, one H
// moved, one fitter asked about both blocks; block 2's forward must be the one a FRESH fitter gives it.  It covers the three molecular
// fitters (the δ scalar fitter, the grid density fitter, the Dunlap fitter); the sampler, the G-space projector and the Hubbard
// projectors take the same InstanceSlotMap change.
#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <string>

import qchem.ChargeDensity;
import qchem.BasisSet;                          // Real_BS, VxcFit, rFIT_SF_ABS, PairProjector
import qchem.BasisSet.Gaussian.Point.Factory;   // Gaussian::Factory
import qchem.BasisSet.CoulombDeltaFit_IBS;      // the δ Coulomb fit basis
import qchem.Fitting.FunctionFitter;            // MakeDeltaScalarFitter
import qchem.Fitting.CoulombOperator;           // MakeDunlapDensityFitter / MakeGridDensityFitter
import qchem.Mesh;                              // qcMesh::MeshParams
import qchem.Mesh.CoulombOperator;              // MakeDirectSumCoulomb
import qchem.Symmetry;                          // SymmetryOrder
import qchem.Symmetry.Unit;                     // UnitQN
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;
import qchem.Structure.MeshFactory;   // MakeMeshMaker: the quadrature a delta fit is built over (V-MESH-FACTORY step 4)
using namespace qchem;

namespace
{
Molecule Water(double dy)
{
    Molecule m;
    m.Insert(Site{8, rvec3_t( 0.0,     0.0,        0.0)});
    m.Insert(Site{1, rvec3_t( 1.43042, 1.10715,    0.0)});
    m.Insert(Site{1, rvec3_t(-1.43042, 1.10715+dy, 0.0)});
    return m;
}
const BasisSet::Real_DFT_OIBS* FirstBlock(const BasisSet::Real_BS& bs)
{
    for (auto b : bs.Iterate<BasisSet::Real_DFT_OIBS>()) return b;
    return nullptr;
}
//! A deterministic PSD density matrix over \a n functions.
hmat_t<double> SomeD(size_t n)
{
    hmat_t<double> D(n);
    for (size_t i=0; i<n; i++) for (size_t j=0; j<=i; j++)
    {
        double s=0.0;
        for (size_t m=0; m<5; m++) s+=std::sin(0.3+1.1*i+0.7*m)*std::sin(0.3+1.1*j+0.7*m);
        D(i,j)=2.0*s/double(n);
    }
    return D;
}
double MaxDiff(const rvec_t& a, const rvec_t& b)
{
    EXPECT_EQ(a.size(), b.size());
    double d=0.0;
    for (size_t g=0; g<std::min(a.size(),b.size()); g++) d=std::max(d, std::fabs(a[g]-b[g]));
    return d;
}
} // namespace

TEST(FitterBlockKeys, TwoBlocksOfEqualSymmetryAskedOfOneFitterGetTheirOwnViews)
{
    using namespace BasisSet::Gaussian;
    const Molecule mol1=Water(0.0), mol2=Water(0.1);   // one H moved 0.1 bohr: the same symmetry, different functions
    const std::unique_ptr<BasisSet::Real_BS> bs1(Factory(BasisSetData::DZVP, &mol1, Engine::MnD, Angular::Cartesian));
    const std::unique_ptr<BasisSet::Real_BS> bs2(Factory(BasisSetData::DZVP, &mol2, Engine::MnD, Angular::Cartesian));
    const auto* b1=FirstBlock(*bs1);
    const auto* b2=FirstBlock(*bs2);
    ASSERT_TRUE(b1 && b2);
    const SymmetryOrder less;
    ASSERT_FALSE(less(b1->GetSymt(), b2->GetSymt()) || less(b2->GetSymt(), b1->GetSymt()))
        << "the two blocks must share a spatial symmetry, or the gate is vacuous";
    const hmat_t<double> D=SomeD(b2->GetNumFunctions());

    const qcMesh::MeshParams mesh{.radial=qcMesh::RadialKind::MHL, .nRadial=40, .mhl_m=3, .mhl_alpha=2.0,
                                  .angular=qcMesh::AngularKind::Lebedev, .angularDegree=17, .beckeOrder=2};
    BasisSet::FitQuadrature quad=MakeMeshMaker(mol1)->Make(mesh, qcMesh::NoSymmetry());
    std::shared_ptr<const BasisSet::FIT_SF_ABS<double>> delta(bs1->CreateVxcFitBasisSet(&mol1, BasisSet::VxcFit::Delta, &quad));
    const auto cdelta=std::make_shared<const BasisSet::CoulombDeltaFit_IBS>(std::make_shared<BasisSet::BasisCache>(), delta,
                                                                             qcMesh::MakeDirectSumCoulomb(*quad.GetMesh()), sym_t(new UnitQN));
    const std::shared_ptr<const BasisSet::rFIT_CD_ABS> gauss(bs1->CreateCDFitBasisSet(&mol1));

    // Each maker builds a fresh fitter over the SAME fit basis; the fitter is held through its forward face.
    struct Case {std::string name; std::function<std::shared_ptr<const BasisSet::PairProjector<double>>()> make;};
    const Case cases[]={
        {"delta scalar", [&]{return std::shared_ptr<const BasisSet::PairProjector<double>>(
                                 std::dynamic_pointer_cast<const BasisSet::PairProjector<double>>(
                                     std::shared_ptr<const Fitting::DeltaScalarFit<double>>(Fitting::MakeDeltaScalarFitter(delta))));}},
        {"grid density", [&]{return std::dynamic_pointer_cast<const BasisSet::PairProjector<double>>(
                                 std::shared_ptr<const Fitting::FirstOrderDensityFit<double>>(Fitting::MakeGridDensityFitter(cdelta, quad.GetMesh())));}},
        {"dunlap",       [&]{return std::dynamic_pointer_cast<const BasisSet::PairProjector<double>>(
                                 std::shared_ptr<const Fitting::DunlapDensityFit<double>>(Fitting::MakeDunlapDensityFitter(gauss)));}},
    };
    for (const Case& c : cases)
    {
        SCOPED_TRACE(c.name);
        const auto fitter=c.make(), fresh=c.make();
        ASSERT_TRUE(fitter && fresh) << "the fitter vends no PairProjector<double>";
        const qcMesh::MatrixForward<double>& f1=fitter->Forward(*b1);
        const qcMesh::MatrixForward<double>& f2=fitter->Forward(*b2);
        EXPECT_NE(&f1, &f2) << "block 2 was served block 1's view";
        const rvec_t p1=f1.Forward(D).Vector(), p2=f2.Forward(D).Vector(), pf=fresh->Forward(*b2).Forward(D).Vector();
        EXPECT_GT(MaxDiff(p2,p1), 1e-6) << "block 2's forward must differ from block 1's (its H moved)";
        EXPECT_EQ(MaxDiff(p2,pf), 0.0) << "block 2's forward must be the one a fresh fitter gives it";
        EXPECT_EQ(&f1, &fitter->Forward(*b1)) << "asked again, block 1 is a HIT";
    }
}
