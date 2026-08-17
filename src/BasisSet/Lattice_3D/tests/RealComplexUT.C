// File: RealComplexUT.C  Unit tests for the Step-3b MIXED basis-set child slot (doc/RealComplexPlan.md).
//
// The whole-set container (BasisSetImp) is the DECISION POINT of the realness rule (§4 item 1): its
// children are typed by their OWN orbital scalar, so a complex-faced periodic set can hold a real TRIM
// block beside general-k complex blocks.  These tests build exactly that mixed set -- a real Γ tGPW_IBS
// <double> and a complex k=(1/4,0,0) tGPW_IBS<dcmplx> in one BasisSetImp<dcmplx> -- and pin the 3b
// contract: the scalar-independent whole-set faces aggregate across BOTH alternatives; the typed child
// view exposes each block with its own scalar; and the face's same-scalar Iterate view THROWS (loud, not
// silent) on the cross-scalar child until Step 3c re-threads the consumers per block.
//
// A unit test may cheat-import Internal modules (CLAUDE.md); production code never subclasses the slot
// this directly.
#include "gtest/gtest.h"
#include <memory>
#include <stdexcept>
#include <variant>

import qchem.BasisSet.Internal.BasisSetImp;   // BasisSetImp<dcmplx> + ibs_child_t (the slot under test)
import qchem.BasisSet.Lattice_3D.GPW_IBS;     // tGPW_IBS<T> (both alternatives)
import qchem.BasisSet.Molecule.Factory;       // the molecular Gaussian basis over the cell's atoms
import qchem.Lattice_3D;                      // UnitCell
import qchem.Types;

using namespace qchem;
using namespace qchem::BasisSet;
using BasisSet::Lattice_3D::tGPW_IBS;

namespace
{
// A minimal concrete over the slot: Insert is protected (concretes own construction), so the test
// subclass republishes it.  Face = Complex_BS, exactly like GPW_BasisSet.
class MixedTestBS : public BasisSetImp<dcmplx>
{
public:
    using BasisSetImp<dcmplx>::Insert;     // both typed overloads
    using BasisSetImp<dcmplx>::GetChild;
    size_t NumIBS() const {return GetNumIBS();}
    const tBasisSet<dcmplx>::obs_t* IBS(size_t i) const {return GetIBS(i);}
};

std::shared_ptr<const Real_BS> MakeMol(const UnitCell& cell)
{
    namespace M = qchem::BasisSet::Molecule;   // disambiguate from qchem::Molecule (the Structure)
    return std::shared_ptr<const Real_BS>(
        M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
}
} //anon

TEST(RealComplexBasisSlot, MixedChildrenAggregateAndTypedViewsWork)
{
    UnitCell cell(8.0);
    cell.AddAtom(14,{0.5,0.5,0.5});
    auto mol = MakeMol(cell);

    MixedTestBS bs;
    auto* re = new tGPW_IBS<double>(cell, ivec3_t(4,4,4), ivec3_t(0,0,0), mol, 0.0);   // Γ: TRIM, REAL
    auto* cx = new tGPW_IBS<dcmplx>(cell, ivec3_t(4,4,4), ivec3_t(1,0,0), mol, 0.0);   // k=(1/4,0,0): complex
    bs.Insert(re);
    bs.Insert(cx);
    ASSERT_EQ(bs.NumIBS(), 2u);

    // Scalar-independent whole-set faces aggregate across BOTH alternatives (generic visits).
    EXPECT_EQ(bs.GetNumFunctions(), re->GetNumFunctions()+cx->GetNumFunctions());
    auto irreps = bs.GetIrreps(Spin::None);
    ASSERT_EQ(irreps.size(), 2u);
    EXPECT_TRUE (irreps[0].sym->IsReal());   // Γ is TRIM (Step 1's exact query)
    EXPECT_FALSE(irreps[1].sym->IsReal());   // (1/4,0,0) is not

    // The typed child view exposes each block with its OWN scalar.
    EXPECT_NE(std::get_if<std::unique_ptr<Orbital_1E_IBS<double>>>(&bs.GetChild(0)), nullptr);
    EXPECT_NE(std::get_if<std::unique_ptr<Orbital_1E_IBS<dcmplx>>>(&bs.GetChild(1)), nullptr);
    EXPECT_TRUE (std::get<std::unique_ptr<Orbital_1E_IBS<double>>>(bs.GetChild(0))->IsReal());

    // The face's SAME-SCALAR view: the complex child is reachable, the real child throws LOUDLY (the
    // 3c contract -- nothing on a shipped path builds a mixed set until the consumers are per-block).
    EXPECT_EQ(bs.IBS(1), static_cast<const tBasisSet<dcmplx>::obs_t*>(cx));
    EXPECT_THROW(bs.IBS(0), std::logic_error);
}

TEST(RealComplexBasisSlot, HomogeneousComplexSetIsUnchanged)
{
    UnitCell cell(8.0);
    cell.AddAtom(14,{0.5,0.5,0.5});
    auto mol = MakeMol(cell);

    MixedTestBS bs;
    auto* g0 = new tGPW_IBS<dcmplx>(cell, ivec3_t(2,2,2), ivec3_t(0,0,0), mol, 0.0);
    auto* g1 = new tGPW_IBS<dcmplx>(cell, ivec3_t(2,2,2), ivec3_t(1,1,1), mol, 0.0);
    bs.Insert(g0);
    bs.Insert(g1);

    // The pre-3b behaviour, bit for bit: every child reachable through the typed face.
    EXPECT_EQ(bs.IBS(0), static_cast<const tBasisSet<dcmplx>::obs_t*>(g0));
    EXPECT_EQ(bs.IBS(1), static_cast<const tBasisSet<dcmplx>::obs_t*>(g1));
    EXPECT_EQ(bs.GetNumFunctions(), g0->GetNumFunctions()+g1->GetNumFunctions());
}
