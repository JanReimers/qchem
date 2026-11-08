// File: src/Structure/tests/ImpositionUT.C  ResolveImposition, the run's symmetry policy (V-BS-SYMMETRY S1): every face is checked
// against its DEFINITION in terms of the lattice's own detection (GetSpaceGroup / ShubnikovOps), over every op -- free, grey
// imposed (Γ and a shifted k mesh), and magnetic (MnO AFM-II).
#include "gtest/gtest.h"
#include <vector>
#include <cmath>

import qchem.Lattice_3D.Imposition;
import qchem.Symmetry.Lattice_3D.Fold;   // MapsMeshOntoItself (the density face's mesh restriction)
import qchem.Types;

using namespace qchem;
namespace SL = qchem::Symmetry::Lattice_3D;

namespace
{
bool Same(const Matrix3D<double>& A, const Matrix3D<double>& B)
{
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (std::fabs(A(i,j)-B(i,j))>1e-12) return false;
    return true;
}
bool Same(const rvec3_t& a, const rvec3_t& b) {return std::fabs(a.x-b.x)+std::fabs(a.y-b.y)+std::fabs(a.z-b.z)<1e-12;}

FCCUnitCell SiDiamond()
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    return cell;
}
//! MnO AFM-II: the FCC cell doubled along [111] (as SiteGroupsUT).
UnitCell MnO_AFM2(double a=8.40)
{
    UnitCell cell=BravaisCell(Bravais::CubicF, {.a=a}, Matrix3D<int>(0,1,1, 1,0,1, 1,1,0));
    cell.AddAtom(25, {0,0,0});
    cell.AddAtom(25, {0.5,0.5,0.5}, true);
    cell.AddAtom(8,  {0.25,0.25,0.25});
    cell.AddAtom(8,  {0.75,0.75,0.75});
    return cell;
}
}

//! A free run imposes nothing; its XC mesh (when asked) is invariant under the WHOLE detected group.
TEST(Imposition, AFreeRunImposesNothingAndKeepsTheDetectedMeshGroup)
{
    const Lattice_3D lat(SiDiamond(), ivec3_t(2,2,2));
    const ImpositionOptions o{.imposeSymmetry=false};
    const CrystalImposition imp=ResolveImposition(lat, o);
    EXPECT_FALSE(imp.policy.Imposes());
    EXPECT_TRUE(imp.options==o);
    EXPECT_TRUE(imp.recipFold.empty() && imp.recipDensity.empty() && imp.directDensity.empty() && imp.magneticDirect.empty());
    const auto& direct=lat.GetSpaceGroup().DirectOps();
    ASSERT_EQ(imp.meshInvariance.size(), direct.size());
    for (size_t i=0;i<direct.size();i++)
    {
        EXPECT_TRUE(Same(imp.meshInvariance[i].W, direct[i].W) && Same(imp.meshInvariance[i].tau, direct[i].tau));
        EXPECT_EQ(imp.meshInvariance[i].sigma, SL::SpinAction::None);
    }
    EXPECT_TRUE(ResolveImposition(lat, {.imposeSymmetry=false, .freeMeshInvariant=false}).meshInvariance.empty())
        << "the plain mesh: no invariance group";
}

//! Grey imposition: the fold is the symmorphic linear group with time reversal; the density faces are the detected {W|τ}
//! restricted to the ops that map the k mesh onto itself -- ALL of them on a Γ-centred mesh, a strict subset on a shifted one.
TEST(Imposition, GreyImpositionIsTheDetectedGroupRestrictedToTheMesh)
{
    for (const rvec3_t shift : {rvec3_t(0,0,0), rvec3_t(0.5,0.5,0.5)})
    {
        const Lattice_3D lat(SiDiamond(), ivec3_t(2,2,2));
        const CrystalImposition imp=ResolveImposition(lat, {.imposeSymmetry=true, .kShift=shift});
        const SL::SpaceGroup& sg=lat.GetSpaceGroup();
        EXPECT_TRUE(imp.policy.Imposes());
        EXPECT_TRUE(imp.magneticDirect.empty() && imp.meshInvariance.empty());
        const auto fold=sg.ReciprocalPointOps(true, true);
        ASSERT_EQ(imp.recipFold.size(), fold.size());
        for (size_t i=0;i<fold.size();i++) EXPECT_TRUE(Same(imp.recipFold[i], fold[i]));
        // The density faces: exactly the detected ops that map the mesh onto itself, in detection order, reciprocal = W^T.
        size_t kept=0;
        for (const auto& op : sg.DirectOps())
        {
            if (!SL::MapsMeshOntoItself(Transpose(op.W), lat.GetLimits(), shift)) continue;
            ASSERT_LT(kept, imp.directDensity.size());
            EXPECT_TRUE(Same(imp.directDensity[kept].W, op.W) && Same(imp.directDensity[kept].tau, op.tau));
            EXPECT_TRUE(Same(imp.recipDensity[kept].U, Transpose(op.W)) && Same(imp.recipDensity[kept].tau, op.tau));
            kept++;
        }
        EXPECT_EQ(imp.directDensity.size(), kept);
        EXPECT_EQ(imp.recipDensity.size(), kept);
        if (shift==rvec3_t(0,0,0)) EXPECT_EQ(kept, sg.DirectOps().size()) << "a Γ-centred mesh respects the whole group";
    }
}

//! Magnetic imposition: the Shubnikov group of the decoration, closed; the density faces are its σ=None subgroup (a sublattice
//! swap is not a symmetry of one channel); no k fold yet (Γ-only bookkeeping).  An all-zero decoration is the grey case.
TEST(Imposition, AStaggeringImposesTheShubnikovGroupAndItsSigmaNoneDensity)
{
    const Lattice_3D lat(MnO_AFM2(), ivec3_t(1,1,1));
    const std::vector<int> afm{+1,-1,0,0};
    const CrystalImposition imp=ResolveImposition(lat, {.imposeSymmetry=true, .siteSpins=afm});
    const auto shub=lat.ShubnikovOps(afm);
    ASSERT_EQ(imp.magneticDirect.size(), shub.size());
    size_t nFlip=0, nNone=0;
    for (size_t i=0;i<shub.size();i++)
    {
        EXPECT_TRUE(Same(imp.magneticDirect[i].W, shub[i].W) && Same(imp.magneticDirect[i].tau, shub[i].tau));
        EXPECT_EQ(imp.magneticDirect[i].sigma, shub[i].sigma);
        if (shub[i].sigma==SL::SpinAction::Flip) {nFlip++; continue;}
        ASSERT_LT(nNone, imp.directDensity.size());
        EXPECT_TRUE(Same(imp.directDensity[nNone].W, shub[i].W) && Same(imp.directDensity[nNone].tau, shub[i].tau));
        EXPECT_TRUE(Same(imp.recipDensity[nNone].U, Transpose(shub[i].W)));
        nNone++;
    }
    EXPECT_GT(nFlip, 0u) << "AFM-II has spin-flip ops: the case this branch exists for";
    EXPECT_EQ(imp.directDensity.size(), nNone);
    EXPECT_EQ(imp.recipDensity.size(), nNone);
    EXPECT_TRUE(imp.recipFold.empty() && imp.meshInvariance.empty());

    const CrystalImposition grey=ResolveImposition(lat, {.imposeSymmetry=true, .siteSpins={0,0,0,0}});
    EXPECT_TRUE(grey.magneticDirect.empty());
    EXPECT_EQ(grey.directDensity.size(), lat.GetSpaceGroup().DirectOps().size()) << "no staggering: the grey group";
}

//! THE INTERSECTION RULE ON THE MAGNETIC BRANCH (R1.0r, user ruling 2026-10-29): the imposed Shubnikov group is intersected with the
//! ops that map the k mesh onto itself, exactly as the spatial branch is -- every face (the Shubnikov set itself, the σ=None density
//! faces) carries only ops the mesh respects, in detection order.  An anisotropic mesh (2x2x1) and a shifted one on the MnO AFM-II
//! cell both drop ops; the Γ-centred isotropic one drops none.
TEST(Imposition, AShubnikovImpositionIsRestrictedToTheMeshToo)
{
    const std::vector<int> afm{+1,-1,0,0};
    struct Mesh { ivec3_t N; rvec3_t shift; bool drops; };
    for (const Mesh m : {Mesh{ivec3_t(2,2,2),rvec3_t(0,0,0),false}, Mesh{ivec3_t(2,2,1),rvec3_t(0,0,0),true},
                         Mesh{ivec3_t(2,2,2),rvec3_t(0.5,0,0),true}})
    {
        const Lattice_3D lat(MnO_AFM2(), m.N);
        const CrystalImposition imp=ResolveImposition(lat, {.imposeSymmetry=true, .siteSpins=afm, .kShift=m.shift});
        const auto shub=lat.ShubnikovOps(afm);
        size_t kept=0, nNone=0, nFlip=0;
        for (const auto& op : shub)
        {
            if (!SL::MapsMeshOntoItself(Transpose(op.W), m.N, m.shift)) continue;
            ASSERT_LT(kept, imp.magneticDirect.size());
            EXPECT_TRUE(Same(imp.magneticDirect[kept].W, op.W) && Same(imp.magneticDirect[kept].tau, op.tau));
            EXPECT_EQ(imp.magneticDirect[kept].sigma, op.sigma);
            kept++;
            if (op.sigma==SL::SpinAction::Flip) {nFlip++; continue;}
            ASSERT_LT(nNone, imp.directDensity.size());
            EXPECT_TRUE(Same(imp.directDensity[nNone].W, op.W) && Same(imp.directDensity[nNone].tau, op.tau));
            EXPECT_TRUE(Same(imp.recipDensity[nNone].U, Transpose(op.W)));
            nNone++;
        }
        EXPECT_EQ(imp.magneticDirect.size(), kept);
        EXPECT_EQ(imp.directDensity.size(), nNone);
        EXPECT_EQ(imp.recipDensity.size(), nNone);
        for (const auto& op : imp.directDensity)
            EXPECT_TRUE(SL::MapsMeshOntoItself(Transpose(op.W), m.N, m.shift)) << "no density op the mesh does not respect";
        if (m.drops) EXPECT_LT(kept, shub.size()) << "this mesh does not respect the whole Shubnikov group";
        else         EXPECT_EQ(kept, shub.size()) << "a Γ-centred isotropic mesh respects the whole group";
        EXPECT_GT(nFlip, 0u) << "the kept set still carries spin-flip ops (the order is still imposed)";
    }
}

//! The imposition BUILDS the run's XC mesh under its own ops (XCMeshSource, V-BS-SYMMETRY S2): the cell's plain mesh when it
//! imposes nothing and asks no invariance (point for point what the cell itself builds); an invariant, orbit-folded mesh when it
//! imposes; and, magnetically, one Shubnikov spin tag per op and a flip-fixed flag per mesh point.
TEST(Imposition, ItBuildsTheXCMeshUnderItsOwnOps)
{
    const qcMesh::MeshParams mp{.nUniform=12, .cellKind=qcMesh::UnitCellKind::Uniform};
    const FCCUnitCell si=SiDiamond();
    const Lattice_3D lat(si, ivec3_t(1,1,1));

    const qcMesh::FoldedMesh plain=ResolveImposition(lat, {.imposeSymmetry=false, .freeMeshInvariant=false}).MakeXCMesh(mp);
    const qcMesh::Mesh ref=si.CreateIntegrationMesh(mp);
    EXPECT_FALSE(plain.HasFold());
    ASSERT_EQ(plain.size(), ref.size());
    for (size_t i=0;i<ref.size();i++)
    {
        EXPECT_TRUE(Same(plain.GetMesh()->Points()[i], ref.Points()[i])) << "point " << i;
        EXPECT_EQ(plain.GetMesh()->Weights()[i], ref.Weights()[i]) << "weight " << i;
    }

    const CrystalImposition grey=ResolveImposition(lat, {.imposeSymmetry=true});
    const qcMesh::FoldedMesh imposed=grey.MakeXCMesh(mp);
    EXPECT_TRUE(imposed.HasFold()) << "an imposed run's density is star-averaged on its mesh";
    EXPECT_LT(imposed.NumOrbits(), imposed.size()) << "the 48 Fd-3m ops fold the mesh";
    EXPECT_FALSE(imposed.IsMagnetic());
    double W=0; const auto& w=imposed.GetMesh()->Weights(); for (size_t i=0;i<w.size();i++) W+=w[i];
    EXPECT_NEAR(W, si.GetCellVolume(), 1e-9*si.GetCellVolume()) << "the invariant mesh still integrates the cell";

    const UnitCell mno=MnO_AFM2();
    const Lattice_3D mlat(mno, ivec3_t(1,1,1));
    const CrystalImposition mag=ResolveImposition(mlat, {.imposeSymmetry=true, .siteSpins={+1,-1,0,0}});
    const qcMesh::FoldedMesh shub=mag.MakeXCMesh(mp);
    EXPECT_TRUE(shub.HasFold() && shub.IsMagnetic());
    EXPECT_EQ(shub.NumSpinOps(), mag.magneticDirect.size()) << "one σ tag per Shubnikov op";
    EXPECT_EQ(shub.NumFlipFixed(), shub.size()) << "one flip-fixed flag per mesh point";
}
