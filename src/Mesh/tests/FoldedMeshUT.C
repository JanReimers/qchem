// File: Mesh/tests/FoldedMeshUT.C  FoldedMesh -- the mesh/fold pairing, and the check that makes it a type.
//
// WHY THIS TYPE EXISTS (user, 2026-09-08).  The pointwise star-average is an EXACT projector *only* because
// the mesh is invariant under the ops the fold was built from.  Before this type the two travelled side by
// side in a fitting struct and every consumer wired one to the other by hand; a fold built for a different
// mesh was undetectable until an index happened to go out of range -- or never, if they all happened to
// land in range, in which case the run just averaged the wrong points together.
//
// So the constructor's CHECK is the feature, and these tests are mostly about the check.
#include "gtest/gtest.h"
#include <memory>
#include <vector>

import qchem.Mesh.Folded;
import qchem.Mesh.Builder;
import qchem.Types;

using namespace qchem;
using qcMesh::FoldedMesh;
namespace SL = qchem::Symmetry::Lattice_3D;

namespace {

//! \a n points with unit weights -- geometry is irrelevant here, only the point COUNT is.
std::shared_ptr<const qcMesh::Mesh> FlatMesh(size_t n)
{
    qcMesh::MeshBuilder b;
    for (size_t i=0;i<n;i++) b.Append(rvec3_t(double(i),0,0), 1.0);
    return std::make_shared<const qcMesh::Mesh>(b.take());
}

//! A fold pairing 2i <-> 2i+1: \a n/2 orbits of two, each represented by its even member.
SL::Fold PairFold(size_t n)
{
    SL::Fold f;
    f.owner.resize(n);
    for (size_t i=0;i<n;i++) f.owner[i]=int(i/2);
    for (size_t o=0;o<n/2;o++)
    {
        f.repRaw.push_back(int(2*o));
        f.starSize.push_back(2);
        f.members.push_back({{int(2*o),0},{int(2*o+1),1}});
    }
    return f;
}

} // anon

//---------------------------------------------------------------------------------------------------
// The default is a FREE RUN and must be a silent no-op: no consumer should have to ask "was symmetry
// imposed?" before calling.  That is what lets the engine call Symmetrize unconditionally.
TEST(FoldedMesh, DefaultIsAFreeRunAndSymmetrizingIsANoOp)
{
    FoldedMesh fm;
    EXPECT_FALSE(fm.HasFold());
    EXPECT_FALSE(fm.IsMagnetic());
    EXPECT_EQ(fm.size(), 0u);
    EXPECT_EQ(fm.GetMesh(), nullptr);

    rvec_t f{1.0, 2.0, 3.0};
    fm.Symmetrize(f);
    EXPECT_DOUBLE_EQ(f[0], 1.0);   // untouched: the projector is the identity
    EXPECT_DOUBLE_EQ(f[1], 2.0);
    EXPECT_DOUBLE_EQ(f[2], 3.0);
}

//---------------------------------------------------------------------------------------------------
TEST(FoldedMesh, StarAverageReplacesEveryOrbitByItsMean)
{
    FoldedMesh fm(FlatMesh(4), PairFold(4));
    ASSERT_TRUE(fm.HasFold());
    EXPECT_EQ(fm.size(), 4u);
    EXPECT_EQ(fm.NumOrbits(), 2u);

    rvec_t f{1.0, 3.0, 10.0, 20.0};
    fm.Symmetrize(f);
    EXPECT_DOUBLE_EQ(f[0], 2.0);    // (1+3)/2
    EXPECT_DOUBLE_EQ(f[1], 2.0);
    EXPECT_DOUBLE_EQ(f[2], 15.0);   // (10+20)/2
    EXPECT_DOUBLE_EQ(f[3], 15.0);

    // IDEMPOTENT -- it is a projector, so a second application must change nothing.
    const rvec_t once=f;
    fm.Symmetrize(f);
    for (size_t i=0;i<f.size();i++) EXPECT_DOUBLE_EQ(f[i], once[i]);
}

//---------------------------------------------------------------------------------------------------
// With no fold, every point is its own star -- so NumOrbits must report the point count, not zero.  It is
// the denominator of the fold FACTOR in the run report, and a zero there would read as an infinite win.
TEST(FoldedMesh, WithoutAFoldEveryPointIsItsOwnOrbit)
{
    FoldedMesh fm(FlatMesh(7), SL::Fold{});
    EXPECT_FALSE(fm.HasFold());
    EXPECT_EQ(fm.size(), 7u);
    EXPECT_EQ(fm.NumOrbits(), 7u);
}

//---------------------------------------------------------------------------------------------------
// ★ THE POINT OF THE TYPE: a fold that does not belong to this mesh is rejected AT CONSTRUCTION, with a
// message, in Release as well as Debug.  Each of these was previously either an assert compiled out under
// NDEBUG or nothing at all.
TEST(FoldedMesh, AFoldBuiltForADifferentMeshIsRejected)
{
    // Wrong point count -- the common case: the mesh was rebuilt (tail-dropped, orbit-filtered) and the
    // fold was not.
    EXPECT_THROW(FoldedMesh(FlatMesh(6), PairFold(4)), std::runtime_error);
    EXPECT_THROW(FoldedMesh(FlatMesh(4), PairFold(6)), std::runtime_error);

    // Right count, but an orbit names a point this mesh does not have.  This is the one that used to be
    // undetectable: the sizes agree, so every size assert passes, and the average silently reads garbage.
    SL::Fold f=PairFold(4);
    f.repRaw[1]=99;
    EXPECT_THROW(FoldedMesh(FlatMesh(4), std::move(f)), std::runtime_error);

    // The odd-field audit must cover every point too.
    EXPECT_THROW(FoldedMesh(FlatMesh(4), PairFold(4), {}, std::vector<char>{1,0}), std::runtime_error);
}

//---------------------------------------------------------------------------------------------------
// A well-formed bundle must NOT throw -- the check has to be a filter, not a wall.
TEST(FoldedMesh, AConsistentBundleIsAccepted)
{
    EXPECT_NO_THROW(FoldedMesh(FlatMesh(4), PairFold(4)));
    EXPECT_NO_THROW(FoldedMesh(FlatMesh(4), PairFold(4), {}, std::vector<char>(4,0)));
    EXPECT_NO_THROW(FoldedMesh(nullptr, SL::Fold{}));           // the free-run default
}

//---------------------------------------------------------------------------------------------------
// The magnetic (Shubnikov) projection: rho is EVEN under the orbit mean, m is ODD under the chi-signed
// one, and m must vanish exactly at the flip-fixed points.  Grey semantics (no spin tags) fall back to
// averaging each channel independently -- which is the S4 negative control, and it ERASES a staggered m.
TEST(FoldedMesh, GreySemanticsAverageEachChannelAndSoEraseAStaggeredMoment)
{
    FoldedMesh grey(FlatMesh(4), PairFold(4));       // no sigmas => grey
    EXPECT_FALSE(grey.IsMagnetic());

    rvec_t rho{1.0, 1.0, 1.0, 1.0}, m{+1.0, -1.0, +1.0, -1.0};
    grey.SymmetrizeSpin(rho, m);
    for (size_t i=0;i<4;i++) EXPECT_DOUBLE_EQ(rho[i], 1.0);
    // the +/- pair sits inside ONE orbit, so the spatial mean cancels it: order destroyed, as designed
    for (size_t i=0;i<4;i++) EXPECT_NEAR(m[i], 0.0, 1e-15);
}
